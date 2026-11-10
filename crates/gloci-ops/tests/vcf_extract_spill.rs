//! `extract_matrix` under a budget so small that the codes spill to the temp
//! directory and the rows are written in blocks (the CLI refuses a budget
//! under 16M, so this is exercised here).
use gloci_io::ReaderOptions;
use gloci_ops::{ExtractMode, ExtractOptions, VcfSelect, extract_matrix, open_selected};
use std::io::Write;
use std::num::NonZeroUsize;

fn matrix(
    vcf: &std::path::Path,
    max_mem: u64,
    tmp: &std::path::Path,
) -> (Vec<u8>, gloci_ops::ExtractSummary) {
    let ro = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Ref,
        invalid: gloci_io::InvalidPolicy::Fail,
        chrom_map: None,
    };
    let select = VcfSelect::default();
    let input = open_selected(&vcf.to_string_lossy(), &ro, &select).unwrap();
    let o = ExtractOptions {
        select,
        mode: ExtractMode::Matrix012,
        skip_invalid: false,
        max_mem,
        tmp_dir: tmp.to_path_buf(),
    };
    let (mut m, mut i, mut p) = (Vec::new(), Vec::new(), Vec::new());
    let sum = extract_matrix(input, &o, &mut m, &mut i, &mut p).unwrap();
    (m, sum)
}

#[test]
fn a_tiny_budget_spills_and_equals_the_in_memory_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let n = 37; // not a multiple of 4: packed rows end in a partial byte
    let names: Vec<String> = (0..n).map(|i| format!("p{i}")).collect();
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        names.join("\t")
    );
    let gts = ["0/0", "0/1", "1/1", "./.", "1", "0|1", "./1", "."];
    for site in 0..400usize {
        text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\t.\t.\tGT", 100 + site));
        for i in 0..n {
            text.push('\t');
            text.push_str(gts[(i * 3 + site * 5 + i * site) % gts.len()]);
        }
        text.push('\n');
    }
    let vcf = dir.path().join("wide.vcf");
    std::fs::File::create(&vcf)
        .unwrap()
        .write_all(text.as_bytes())
        .unwrap();
    let spill = dir.path().join("spill");
    std::fs::create_dir(&spill).unwrap();
    let (whole, sum) = matrix(&vcf, 1 << 30, &spill);
    assert!(!sum.spilled && sum.row_blocks == 1, "{sum:?}");
    assert_eq!(std::fs::read_dir(&spill).unwrap().count(), 0);
    // 1 KiB: the codes exceed half of it, so they spill, and a block holds at
    // most one sample's row.
    for budget in [1 << 10, 6 << 10] {
        let (got, sum) = matrix(&vcf, budget, &spill);
        assert_eq!(got, whole, "budget {budget}");
        // it spilled, and the rows were written in several blocks
        assert!(
            sum.spilled && sum.row_blocks > 1,
            "budget {budget}: {sum:?}"
        );
        // the spill file is gone
        assert_eq!(std::fs::read_dir(&spill).unwrap().count(), 0);
    }
}
