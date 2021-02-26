from os import listdir

url = ''
files = listdir()
names = list(set(["_".join(f.split("_")[:-1]) for f in files if "fasta" in f]))
for n in names:
    print(n)
    for t in ["actual", "RNAfold", "rollout"]:
        url = 'http://rna.tbi.univie.ac.at/forna/forna.html?id=fasta&file='

        with open(n+"_"+t+".fasta", 'r') as f:
            text = f.read().replace(' ', '\n')
            url += repr(text)[1:-3]

        if t != "actual":
            url += "&colors="
            with open(n+"_"+t+"_colormap.col", 'r') as f:
                text = f.read()
            
                url += repr(text)[1:-1]
        print(url)
        print('')
    print('')