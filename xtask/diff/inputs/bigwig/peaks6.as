table peaks6
"BED6 plus a signal and a tag"
(
    string chrom; "Chromosome"
    uint chromStart; "Start, 0-based"
    uint chromEnd; "End, exclusive"
    string name; "Name"
    uint score; "Score"
    char[1] strand; "Strand"
    float signal; "Signal value"
    string tag; "Free text tag"
)
