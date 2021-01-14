import hashlib

p1 = "./"
p2 = "../STD/"
files = ["preStatsKmer", "postStatsKmer", "mDuplicateCount", "mDuplicateGC", "mDuplicateDups", "p.fq"]

for it in files:

    try:
        f1 = p1 + it
        with open(f1, 'rb') as fp:
            data = fp.read()
        m1 = hashlib.md5(data).hexdigest()
        # print(m1)
        f2 = p2 + it
        with open(f2, 'rb') as fp:
            data = fp.read()
        m2 = hashlib.md5(data).hexdigest()
        # print(m2)
        if m1 != m2:
            print("GG on check %s" % it)
            ok = 0
        else:
            print("pass on test %s" % it)
    except:
        print("some errors occur when open %s" % it)
        continue
