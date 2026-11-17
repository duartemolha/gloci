//! `vcf query`'s per-sample evaluation (D172): each node gives a site value or
//! a pass per selected sample, combined as bcftools query combines them.

use super::eval::{Cx, need_gt, node_holds, sample_passes};
use super::{EvalScratch, Field, Node, QueryPass};
use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::SiteStats;

/// The value of a node: one truth value for the record, or a pass per header
/// sample (false for unselected samples) with whether the node holds. A node
/// holds when a sample passes, or when it holds through a site value over a
/// record with no selected sample (so `QUAL>8 || GT="het"` keeps a sites-only
/// record, as bcftools does).
enum Q {
    Site(bool),
    Mask(Vec<bool>, bool),
}

/// Does the node give per-sample passes? Static: a `GT`, `FMT/DP` or `FMT/GQ`
/// comparison outside `N_PASS`, or a list or negation holding one.
fn per_sample(n: &Node) -> bool {
    match n {
        Node::And(v) | Node::Or(v) => v.iter().any(per_sample),
        Node::Not(inner) => per_sample(inner),
        Node::Cmp { field, .. } => {
            matches!(field, Field::Gt(_) | Field::FmtDp(_) | Field::FmtGq(_))
        }
        Node::Truthy(_) => false,
    }
}

pub(super) fn eval_query(
    root: &Node,
    exclude: bool,
    site: &SiteView,
    gt: Option<&dyn GenotypeSource>,
    stats: Option<&SiteStats>,
    scratch: &mut EvalScratch,
    pass: &mut Vec<bool>,
) -> Result<QueryPass, String> {
    let mut cx = Cx {
        site,
        gt,
        stats,
        scratch,
        sample: None,
    };
    Ok(match node(root, &mut cx)? {
        Q::Site(holds) => {
            if holds != exclude {
                QueryPass::All
            } else {
                QueryPass::Drop
            }
        }
        Q::Mask(mut m, holds) => {
            let out = if !exclude {
                if holds {
                    QueryPass::Samples
                } else {
                    QueryPass::Drop
                }
            } else if holds {
                flip(&mut m, &cx);
                if m.iter().any(|&b| b) {
                    QueryPass::Samples
                } else {
                    QueryPass::Drop
                }
            } else {
                QueryPass::All
            };
            if out == QueryPass::Samples {
                pass.clear();
                pass.extend_from_slice(&m);
            }
            cx.scratch.masks.push(m);
            out
        }
    })
}

fn selected(cx: &Cx, i: usize) -> bool {
    cx.scratch
        .keep
        .as_ref()
        .is_none_or(|k| k.get(i).copied().unwrap_or(true))
}

/// Flip each selected sample's pass; unselected ones stay false.
fn flip(m: &mut [bool], cx: &Cx) {
    for (i, b) in m.iter_mut().enumerate() {
        *b = selected(cx, i) && !*b;
    }
}

fn set_all_selected(m: &mut [bool], cx: &Cx) {
    for (i, b) in m.iter_mut().enumerate() {
        *b = selected(cx, i);
    }
}

/// A mask of `n` samples, all false, from the pool.
fn blank(cx: &mut Cx, n: usize) -> Vec<bool> {
    let mut m = cx.scratch.masks.pop().unwrap_or_default();
    m.clear();
    m.resize(n, false);
    m
}

fn node(n: &Node, cx: &mut Cx) -> Result<Q, String> {
    if !per_sample(n) {
        // A site-only subtree is D71's rule whole: same short-circuit, same errors.
        return Ok(Q::Site(node_holds(n, cx)?));
    }
    match n {
        Node::And(v) | Node::Or(v) => {
            let is_and = matches!(n, Node::And(_));
            let mut acc = node(&v[0], cx)?;
            for next in &v[1..] {
                // Once the result cannot change, the rest is not read.
                let done = match &acc {
                    Q::Site(b) => *b != is_and,
                    Q::Mask(m, holds) => {
                        if is_and {
                            !holds
                        } else {
                            *holds && m.iter().enumerate().all(|(i, &b)| b || !selected(cx, i))
                        }
                    }
                };
                if done {
                    break;
                }
                let r = node(next, cx)?;
                acc = if is_and {
                    and(acc, r, cx)
                } else {
                    or(acc, r, cx)
                };
            }
            // A mixed list is per sample whichever operand ended it.
            if let Q::Site(b) = acc {
                let count = need_gt(cx)?.n_samples();
                let mut m = blank(cx, count);
                if b {
                    set_all_selected(&mut m, cx);
                }
                acc = Q::Mask(m, b);
            }
            Ok(acc)
        }
        Node::Not(inner) => Ok(match node(inner, cx)? {
            Q::Site(b) => Q::Site(!b),
            Q::Mask(mut m, _) => {
                flip(&mut m, cx);
                let holds = m.iter().any(|&b| b);
                Q::Mask(m, holds)
            }
        }),
        Node::Cmp {
            field: field @ (Field::Gt(sub) | Field::FmtDp(sub) | Field::FmtGq(sub)),
            op,
            rhs,
        } => {
            let src = need_gt(cx)?;
            let count = src.n_samples();
            let mut m = blank(cx, count);
            let range = match *sub {
                Some(k) if k >= count => {
                    return Err(format!(
                        "sample index {k} is too large: the record has {count} samples"
                    ));
                }
                Some(k) => k..k + 1,
                None => 0..count,
            };
            for i in range {
                if !selected(cx, i) {
                    continue;
                }
                src.sample(i, &mut cx.scratch.fields)?;
                let s = &mut *cx.scratch;
                m[i] = sample_passes(field, *op, rhs, &mut s.gt_text, &s.fields);
            }
            let holds = m.iter().any(|&b| b);
            Ok(Q::Mask(m, holds))
        }
        _ => unreachable!("only list, negation and sample comparisons are per sample"),
    }
}

fn and(a: Q, b: Q, cx: &mut Cx) -> Q {
    match (a, b) {
        (Q::Site(a), Q::Site(b)) => Q::Site(a && b),
        (Q::Site(h), Q::Mask(mut m, mh)) | (Q::Mask(mut m, mh), Q::Site(h)) => {
            if !h {
                m.fill(false);
            }
            Q::Mask(m, h && mh)
        }
        (Q::Mask(mut a, ah), Q::Mask(b, bh)) => {
            let both = ah && bh;
            if both {
                for (x, y) in a.iter_mut().zip(&b) {
                    *x |= *y;
                }
            } else {
                a.fill(false);
            }
            cx.scratch.masks.push(b);
            Q::Mask(a, both)
        }
    }
}

fn or(a: Q, b: Q, cx: &mut Cx) -> Q {
    match (a, b) {
        (Q::Site(a), Q::Site(b)) => Q::Site(a || b),
        (Q::Site(h), Q::Mask(mut m, mh)) | (Q::Mask(mut m, mh), Q::Site(h)) => {
            if h {
                set_all_selected(&mut m, cx);
            }
            Q::Mask(m, h || mh)
        }
        (Q::Mask(mut a, ah), Q::Mask(b, bh)) => {
            let either = ah || bh;
            if either {
                set_all_selected(&mut a, cx);
            } else {
                a.fill(false);
            }
            cx.scratch.masks.push(b);
            Q::Mask(a, either)
        }
    }
}
