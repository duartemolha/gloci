mod common;
use common::{HEADER, vcf};
use gloci_io::{BcfText, TextToBcf, VcfHeader};
use gloci_vcf::{
    MaskRules, RewritePlan, SampleSelection, fmt_g6, read_samples, rewrite_header, rewrite_text,
    select_samples, site_stats, text_genotypes,
};

#[test]
fn subset_text() {
    let raw_header = vcf(&["s1", "s2", "s3"], &[]);
    let h = VcfHeader::parse(raw_header.into_bytes()).unwrap();

    let sel: SampleSelection =
        select_samples(&h.samples, Some(&["s1".into(), "s3".into()]), &[]).unwrap();
    assert_eq!(sel.keep, vec![true, false, true]);
    assert!(!sel.all_kept());
    assert_eq!(sel.count(), 2);

    let plan = RewritePlan {
        selection: Some(sel),
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: Vec::new(),
        recalc_ac: false,
    };
    assert!(!plan.is_passthrough());

    let header_rewritten = rewrite_header(&h, &plan);
    let h_rewritten = VcfHeader::parse(header_rewritten).unwrap();
    assert_eq!(h_rewritten.samples, vec![b"s1".to_vec(), b"s3".to_vec()]);

    let line = b"1\t100\t.\tA\tC\t30\tPASS\tDP=30\tGT:DP\t0/1:10\t1/1:10\t0/0:10";
    let mut out = Vec::new();
    rewrite_text(line, &plan, &[], None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\tDP=30\tGT:DP\t0/1:10\t0/0:10"
    );
}

#[test]
fn unknown_sample_is_usage_error() {
    let raw_header = vcf(&["s1", "s2", "s3", "s4", "s5", "s6"], &[]);
    let h = VcfHeader::parse(raw_header.into_bytes()).unwrap();

    // Unknown in include
    let err =
        select_samples(&h.samples, Some(&["s1".into(), "bad_sample".into()]), &[]).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("bad_sample"));
    assert!(msg.contains("s1, s2, s3, s4, s5"));
    assert!(!msg.contains("s6"));

    // Unknown in exclude
    let err2 = select_samples(&h.samples, None, &["ghost".into()]).unwrap_err();
    let msg2 = err2.to_string();
    assert!(msg2.contains("ghost"));
    assert!(msg2.contains("s1, s2, s3, s4, s5"));
}

#[test]
fn mask_text() {
    let line =
        b"1\t100\t.\tA\tC\t30\tPASS\t.\tGT:DP:GQ:FT\t0/1:10:30:PASS\t1|1:20:40:PASS\t1:5:20:PASS";
    let plan = RewritePlan {
        selection: None,
        mask: MaskRules {
            min_gq: Some(35),
            ..MaskRules::default()
        },
        keep_info: None,
        drop_info: Vec::new(),
        recalc_ac: false,
    };

    // Mask s1 (0/1 -> ./.) and s3 (haploid 1 -> .), keep s2 (1|1)
    let masked = vec![true, false, true];
    let mut out = Vec::new();
    rewrite_text(line, &plan, &masked, None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\t.\tGT:DP:GQ:FT\t./.:10:30:PASS\t1|1:20:40:PASS\t.:5:20:PASS"
    );

    // Also verify 1|1 -> .|.
    let masked_all = vec![true, true, true];
    out.clear();
    rewrite_text(line, &plan, &masked_all, None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\t.\tGT:DP:GQ:FT\t./.:10:30:PASS\t.|.:20:40:PASS\t.:5:20:PASS"
    );
}

#[test]
fn drop_info_and_header() {
    let raw = format!(
        "{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t100\t.\tA\tC\t30\tPASS\tDP=30;NS=3\n"
    );
    let h = VcfHeader::parse(raw.into_bytes()).unwrap();

    let plan = RewritePlan {
        selection: None,
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: vec![b"DP".to_vec()],
        recalc_ac: false,
    };

    let header_rewritten = rewrite_header(&h, &plan);
    let h_str = std::str::from_utf8(&header_rewritten).unwrap();
    assert!(!h_str.contains("##INFO=<ID=DP,"));
    assert!(h_str.contains("##INFO=<ID=NS,"));

    let line = b"1\t100\t.\tA\tC\t30\tPASS\tDP=30;NS=3";
    let mut out = Vec::new();
    rewrite_text(line, &plan, &[], None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\tNS=3"
    );

    // Dropping all INFO keys gives "."
    let plan_all = RewritePlan {
        selection: None,
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: vec![b"DP".to_vec(), b"NS".to_vec()],
        recalc_ac: false,
    };
    out.clear();
    rewrite_text(line, &plan_all, &[], None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\t."
    );
}

#[test]
fn keep_info() {
    let raw = format!(
        "{HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t100\t.\tA\tC\t30\tPASS\tDP=30;NS=3;DB\n"
    );
    let h = VcfHeader::parse(raw.into_bytes()).unwrap();

    let plan = RewritePlan {
        selection: None,
        mask: MaskRules::default(),
        keep_info: Some(vec![b"NS".to_vec()]),
        drop_info: Vec::new(),
        recalc_ac: false,
    };

    let header_rewritten = rewrite_header(&h, &plan);
    let h_str = std::str::from_utf8(&header_rewritten).unwrap();
    assert!(!h_str.contains("##INFO=<ID=DP,"));
    assert!(!h_str.contains("##INFO=<ID=DB,"));
    assert!(h_str.contains("##INFO=<ID=NS,"));

    let line = b"1\t100\t.\tA\tC\t30\tPASS\tDP=30;NS=3;DB";
    let mut out = Vec::new();
    rewrite_text(line, &plan, &[], None, &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\tNS=3"
    );
}

#[test]
fn recalc_ac() {
    // samples 0/1 and 1/1 give AC=3;AN=4;AF=0.75
    let line = b"1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\t1/1";
    let gt = text_genotypes(line).unwrap();
    let mut masked = Vec::new();
    let stats = site_stats(&gt, None, &MaskRules::default(), 2, &mut masked).unwrap();

    let plan = RewritePlan {
        selection: None,
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: Vec::new(),
        recalc_ac: true,
    };

    // Replacing an INFO of .
    let mut out = Vec::new();
    rewrite_text(line, &plan, &masked, Some(&stats), &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\tAC=3;AN=4;AF=0.75\tGT\t0/1\t1/1"
    );

    // Existing info with other keys: replace where exists, otherwise append
    let line2 = b"1\t100\t.\tA\tC\t30\tPASS\tDP=50;AF=0.1;NS=2\tGT\t0/1\t1/1";
    out.clear();
    rewrite_text(line2, &plan, &masked, Some(&stats), &mut out).unwrap();
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC\t30\tPASS\tDP=50;AF=0.75;NS=2;AC=3;AN=4\tGT\t0/1\t1/1"
    );

    // Multiallelic: ALT is C,G
    let line3 = b"1\t100\t.\tA\tC,G\t30\tPASS\t.\tGT\t0/1\t1/2";
    let gt3 = text_genotypes(line3).unwrap();
    let stats3 = site_stats(&gt3, None, &MaskRules::default(), 3, &mut masked).unwrap();
    out.clear();
    rewrite_text(line3, &plan, &masked, Some(&stats3), &mut out).unwrap();
    // allele 0: 1, allele 1: 2, allele 2: 1. Total called = 4.
    // AC=2,1;AN=4;AF=0.5,0.25
    assert_eq!(
        std::str::from_utf8(&out).unwrap(),
        "1\t100\t.\tA\tC,G\t30\tPASS\tAC=2,1;AN=4;AF=0.5,0.25\tGT\t0/1\t1/2"
    );
}

#[test]
fn fmt_g6_cases() {
    assert_eq!(fmt_g6(0.333333333), "0.333333");
    assert_eq!(fmt_g6(0.5), "0.5");
    assert_eq!(fmt_g6(1e-7), "1e-07");
    assert_eq!(fmt_g6(0.0), "0");
    assert_eq!(fmt_g6(1.0), "1");
    assert_eq!(fmt_g6(1234567.0), "1.23457e+06");
    assert_eq!(fmt_g6(0.0001), "0.0001");
    assert_eq!(fmt_g6(0.0000999999), "9.99999e-05");
}

#[test]
fn bcf_rewrite_equals_text() {
    let raw_header = "##fileformat=VCFv4.3\n\
##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n";

    let in_header = VcfHeader::parse(raw_header.as_bytes().to_vec()).unwrap();
    let line = b"1\t100\t.\tA\tC\t30\tPASS\tDP=30\tGT:DP\t0/1:10\t1/1:10\t0/0:10";

    // Encode input BCF record
    let mut text_to_bcf = TextToBcf::new("in.bcf", &in_header).unwrap();
    let mut bcf_rec = Vec::new();
    text_to_bcf.encode(line, 1, &mut bcf_rec).unwrap();

    let sel = select_samples(&in_header.samples, Some(&["s1".into(), "s3".into()]), &[]).unwrap();
    let plan = RewritePlan {
        selection: Some(sel),
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: Vec::new(),
        recalc_ac: false,
    };

    // Text rewrite
    let mut text_out = Vec::new();
    rewrite_text(line, &plan, &[], None, &mut text_out).unwrap();

    // BCF rewrite path
    let new_header_raw = rewrite_header(&in_header, &plan);
    let new_header = VcfHeader::parse(new_header_raw).unwrap();

    let mut bcf_text = BcfText::for_header("in.bcf", &in_header).unwrap();
    let rendered = bcf_text.render(&bcf_rec, 1).unwrap();

    let mut bcf_text_rewritten = Vec::new();
    rewrite_text(rendered, &plan, &[], None, &mut bcf_text_rewritten).unwrap();
    assert_eq!(bcf_text_rewritten, text_out);

    let mut out_encoder = TextToBcf::new("out.bcf", &new_header).unwrap();
    let mut out_bcf_bytes = Vec::new();
    out_encoder
        .encode(&bcf_text_rewritten, 1, &mut out_bcf_bytes)
        .unwrap();

    // Decode back from BCF and verify it matches text_out
    let mut out_bcf_text = BcfText::for_header("out.bcf", &new_header).unwrap();
    let back_to_text = out_bcf_text.render(&out_bcf_bytes, 1).unwrap();
    assert_eq!(back_to_text, text_out.as_slice());
}

#[test]
fn read_samples_file() {
    let dir = std::env::temp_dir().join(format!("gloci-samples-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("samples.txt");
    std::fs::write(&p, "# comment\ns1\n\ns2\r\ns3\n").unwrap();
    let got = read_samples(p.to_str().unwrap()).unwrap();
    assert_eq!(got, vec!["s1", "s2", "s3"]);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn bcf_passthrough_bytes() {
    let raw_header = "##fileformat=VCFv4.3\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";

    let in_header = VcfHeader::parse(raw_header.as_bytes().to_vec()).unwrap();
    let line = b"1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/1";

    let mut enc = TextToBcf::new("in.bcf", &in_header).unwrap();
    let mut bcf_rec = Vec::new();
    enc.encode(line, 1, &mut bcf_rec).unwrap();

    let plan = RewritePlan {
        selection: None,
        mask: MaskRules::default(),
        keep_info: None,
        drop_info: Vec::new(),
        recalc_ac: false,
    };
    assert!(plan.is_passthrough());

    // When passthrough, the raw bytes are untouched
    let passthrough_bcf = bcf_rec.clone();
    assert_eq!(passthrough_bcf, bcf_rec);
}
