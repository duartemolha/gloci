//! `vcf query`'s per-sample evaluation (D172): each node gives a site value or
//! a pass per selected sample, combined as bcftools query combines them.

use super::eval::{Cx, need_gt, node_holds, sample_passes};
use super::{EvalScratch, Field, Node, QueryPass};
use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::SiteStats;

/// The value of a node: one truth value for the record, or a pass per header
/// sample (false for unselected samples).
enum Q {
    Site(bool),
    Mask(Vec<bool>),
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
        Q::Mask(mut m) => {
            let passes = m.iter().any(|&b| b);
            let out = if !exclude {
                if passes {
                    QueryPass::Samples
                } else {
                    QueryPass::Drop
                }
            } else if passes {
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

fn node(n: &Node, cx: &mut Cx) -> Result<Q, String> {
    match n {
        Node::And(v) | Node::Or(v) => {
            let is_and = matches!(n, Node::And(_));
            let mut it = v.iter();
            let mut acc = node(it.next().expect("a list has operands"), cx)?;
            for next in it {
                let r = node(next, cx)?;
                acc = if is_and {
                    and(acc, r, cx)
                } else {
                    or(acc, r, cx)
                };
            }
            Ok(acc)
        }
        Node::Not(inner) => Ok(match node(inner, cx)? {
            Q::Site(b) => Q::Site(!b),
            Q::Mask(mut m) => {
                flip(&mut m, cx);
                Q::Mask(m)
            }
        }),
        Node::Cmp {
            field: field @ (Field::Gt(sub) | Field::FmtDp(sub) | Field::FmtGq(sub)),
            op,
            rhs,
        } => {
            let src = need_gt(cx)?;
            let count = src.n_samples();
            let mut m = cx.scratch.masks.pop().unwrap_or_default();
            m.clear();
            m.resize(count, false);
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
            Ok(Q::Mask(m))
        }
        _ => Ok(Q::Site(node_holds(n, cx)?)),
    }
}

fn and(a: Q, b: Q, cx: &mut Cx) -> Q {
    match (a, b) {
        (Q::Site(a), Q::Site(b)) => Q::Site(a && b),
        (Q::Site(holds), Q::Mask(mut m)) | (Q::Mask(mut m), Q::Site(holds)) => {
            if !holds {
                m.fill(false);
            }
            Q::Mask(m)
        }
        (Q::Mask(mut a), Q::Mask(b)) => {
            if a.iter().any(|&x| x) && b.iter().any(|&x| x) {
                for (x, y) in a.iter_mut().zip(&b) {
                    *x |= *y;
                }
            } else {
                a.fill(false);
            }
            cx.scratch.masks.push(b);
            Q::Mask(a)
        }
    }
}

fn or(a: Q, b: Q, cx: &mut Cx) -> Q {
    match (a, b) {
        (Q::Site(a), Q::Site(b)) => Q::Site(a || b),
        (Q::Site(holds), Q::Mask(mut m)) | (Q::Mask(mut m), Q::Site(holds)) => {
            if holds {
                set_all_selected(&mut m, cx);
            }
            Q::Mask(m)
        }
        (Q::Mask(mut a), Q::Mask(b)) => {
            if a.iter().any(|&x| x) || b.iter().any(|&x| x) {
                set_all_selected(&mut a, cx);
            } else {
                a.fill(false);
            }
            cx.scratch.masks.push(b);
            Q::Mask(a)
        }
    }
}
