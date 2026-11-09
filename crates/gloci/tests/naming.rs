mod common;
use common::*;
use predicates::str::contains;

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }

    fn file(&self, name: &str, text: &str) -> String {
        plain_file(self.dir.path(), name, text)
            .to_string_lossy()
            .into_owned()
    }

    fn bgzf_indexed(&self, name: &str, text: &str) -> String {
        let p = bgzf_file(self.dir.path(), name, text);
        gloci()
            .args(["index", p.to_str().unwrap()])
            .assert()
            .success();
        p.to_string_lossy().into_owned()
    }
}

#[test]
fn chr_prefix_mismatch_warns() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}

#[test]
fn disjoint_names_without_prefix() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "chr2\t15\t25\n");
    let assert = gloci()
        .args(["intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: chr2)"
        )));
    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(!stderr.contains("the names differ by a \"chr\" prefix"));
}

#[test]
fn g_alias_hint() {
    let fx = Fx::new();
    let x = fx.file("x.bed", "chr1\t10\t20\n");
    gloci()
        .args(["sort", &x, "-g", "GRCh38"])
        .assert()
        .code(1)
        .stderr(contains("did you mean -g hg38?"));

    let y = fx.file("y.bed", "1\t10\t20\n");
    gloci()
        .args(["sort", &y, "-g", "hg38"])
        .assert()
        .code(1)
        .stderr(contains("did you mean -g GRCh38?"));

    // Chromosome not present in either genome does not suggest alias
    let z = fx.file("z.bed", "chrZ\t10\t20\n");
    let assert = gloci()
        .args(["sort", &z, "-g", "GRCh38"])
        .assert()
        .code(1)
        .stderr(contains("is not in the built-in genome GRCh38"));
    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    assert!(!stderr.contains("did you mean -g"));
}

#[test]
fn no_warning_when_shared() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "chr1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout("chr1\t10\t20\n")
        .stderr("");
}

#[test]
fn merge_warns_on_disjoint_chroms() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "1\t15\t25\n");
    gloci()
        .args(["merge", &a, &b])
        .assert()
        .success()
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}

#[test]
fn complement_warns_on_disjoint_chroms() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "1\t15\t25\n");
    gloci()
        .args(["complement", &a, &b])
        .assert()
        .success()
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}

#[test]
fn indexed_relation_warns_on_disjoint_chroms() {
    let fx = Fx::new();
    let a = fx.bgzf_indexed("a.bed.gz", "chr1\t10\t20\n");
    let b = fx.bgzf_indexed("b.bed.gz", "1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}

#[test]
fn relation_sort_names_user_files() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t20\t30\nchr1\t10\t15\n");
    let b = fx.file("b.bed", "1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b, "--sort"])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}

#[test]
fn relation_unsorted_warns() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\n");
    let b = fx.file("b.bed", "1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b, "--unsorted"])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "warning: {a} and {b} share no chromosome names ({a}: chr1; {b}: 1); \
the names differ by a \"chr\" prefix — rename one file's chromosomes so they match"
        )));
}
