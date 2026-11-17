import pyBigWig
bw = pyBigWig.open("pybigwig.bw", "w")
bw.addHeader([("chr1", 5000), ("chr2", 3000), ("chr10", 2000), ("chrM", 100)])
bw.addEntries(["chr1", "chr1", "chr1"], [0, 100, 250], ends=[100, 200, 300], values=[1.5, 0.25, -2.0])
bw.addEntries("chr2", [10, 20, 30], values=[3.0, 4.5, 0.1], span=5)
bw.addEntries("chr10", 100, values=[7.0, 8.0, 9.0], span=10, step=20)
bw.close()
