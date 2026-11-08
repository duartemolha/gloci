//! `ITree64` (D40): the unsorted-mode index for chromosomes whose
//! coordinates do not fit coitrees' `i32`.

use gloci_core::Interval;

/// A static interval index over `u64` coordinates: an implicit augmented
/// interval tree, like cgranges. The items are sorted by start in one
/// array, with no pointers. The root of each range of items is its middle
/// item, which also holds the largest end in the range; ranges of at most
/// [`LEAF`] items are scanned instead.
///
/// It returns **candidates**: an item `[start, end)` counts as the closed
/// range `[start, end]`, one base wider (so a point `[p, p)` is `[p, p]`,
/// never empty), and a query returns every item whose closed range meets
/// the query's. The caller re-tests them with the exact rule (spec §4.3).
#[derive(Debug, Clone, Default)]
pub struct ITree64 {
    nodes: Vec<Node>,
}

#[derive(Debug, Clone, Copy)]
struct Node {
    start: u64,
    end: u64,
    /// For a range root, the largest end in its range.
    max_end: u64,
    id: u32,
}

/// Ranges this small are scanned, not split.
const LEAF: usize = 16;

impl ITree64 {
    /// Index `items`, each an interval and the id a query returns for it.
    pub fn new(items: Vec<(Interval, u32)>) -> ITree64 {
        let mut nodes: Vec<Node> = items
            .into_iter()
            .map(|(iv, id)| Node {
                start: iv.start,
                end: iv.end,
                max_end: iv.end,
                id,
            })
            .collect();
        nodes.sort_unstable_by_key(|n| (n.start, n.id));
        augment(&mut nodes);
        ITree64 { nodes }
    }

    /// Append to `out` the id of every item whose closed range
    /// `[start, end]` meets the closed range `[lo, hi]` (`start <= hi` and
    /// `end >= lo`), each once, in no particular order.
    pub fn query(&self, lo: u64, hi: u64, out: &mut Vec<u32>) {
        query(&self.nodes, lo, hi, out);
    }
}

/// Set the root of every range of `nodes` (recursively) to the largest end
/// in its range, and return the largest end in `nodes`.
fn augment(nodes: &mut [Node]) -> u64 {
    if nodes.len() <= LEAF {
        return nodes.iter().map(|n| n.end).max().unwrap_or(0);
    }
    let mid = nodes.len() / 2;
    let (left, rest) = nodes.split_at_mut(mid);
    let (root, right) = rest
        .split_first_mut()
        .expect("a range longer than LEAF has a middle");
    root.max_end = augment(left).max(root.end).max(augment(right));
    root.max_end
}

fn query(nodes: &[Node], lo: u64, hi: u64, out: &mut Vec<u32>) {
    if nodes.len() <= LEAF {
        for n in nodes {
            if n.start > hi {
                break;
            }
            if n.end >= lo {
                out.push(n.id);
            }
        }
        return;
    }
    let mid = nodes.len() / 2;
    let root = &nodes[mid];
    if root.max_end < lo {
        // Everything in this range ends before the query.
        return;
    }
    query(&nodes[..mid], lo, hi, out);
    if root.start > hi {
        // So does every item after it.
        return;
    }
    if root.end >= lo {
        out.push(root.id);
    }
    query(&nodes[mid + 1..], lo, hi, out);
}
