#coding:utf-8

"""Population-specific genome generator using 1000Genome database.
"""

import os, sys, re, subprocess, gzip, collections, random, argparse

complementary = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
gmargin = 1000
vmargin = 100000

# artificial phenotypes
preset_for_selected_phenotypes_by_endo = {
    'rs1042602':'A',#	A;A,freckless	C;C/C;A,normal
    'rs26722':'C', #(TYR)	T;T/T;C,dark	C;C,normal
    'rs16891982':'G', #	G;G,light skin	G;C/C;C,black hair
    'rs1426654':'A',#	A;A,light skin	A;G/G;G,dark
    'rs642742':'G', #	G,ligher	A,darker
    'rs12821256':'C',#	C;C,blond
    # music talent
    'rs3057':None,
    'rs2028030':None,
    'rs1007750':None,
    'rs2169325':None,
    'rs1482308':None,
    'rs9854612':None,
    'rs13146789':None,
    'rs12510781':None,
}
preset_snps = {} # set this variable as 'preset_for_selected_phenotypes_by_endo' to use some phenotype-assoociated SNPs

basedir = '.'
genome = os.path.join(basedir, 'hs37d5.fa.gz')
descfile = os.path.join(basedir, 'integrated_call_samples_v3.20130502.ALL.panel')

parser = argparse.ArgumentParser()
parser.add_argument('-g', default=genome, metavar='filename', help='gzipped genome')
parser.add_argument('-s', default=descfile, help='SNP panel', metavar='filename')
parser.add_argument('-p', metavar='filename', default='JPT', help='population group')
parser.add_argument('-o', default='jpn.fa.gz', metavar='filename', help='gzipped output filename')
parser.add_argument('-d', default=basedir, metavar='directory', help='VCF file container directory')
parser.add_argument('--lowercase', action='store_true', help='replaced nucleotides are shown in lower case')
parser.add_argument('--male', action='store_true', help='output genome is set as male (including chrY)')
parser.add_argument('--verbose', action='store_true', help='verbosity')

args = parser.parse_args()

if args.verbose:
    sys.stderr.write('Genome file (gzipped)      : {}\n'.format(args.g))
    sys.stderr.write('SNP panel (.panel file)    : {}\n'.format(args.s))
    sys.stderr.write('Gzipped VCF file directory : {}\n'.format(args.d))
    sys.stderr.write('Population                 : {}\n'.format(args.p))
    sys.stderr.write('Gzipped output file        : {}\n'.format(args.o))

varfiles = {}
female_columns = collections.OrderedDict()
basedir = args.d
genome = args.g
descfile = args.s

# variation files for chromosomes
for fn in os.listdir(basedir):
    m = re.match('ALL.(chr\\w+).phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', fn)
    if m: varfiles[m.group(1)] = os.path.join(basedir, fn)

# population
with open(descfile) as fi:
    fi.readline()
    for line in fi:
        items = line.strip().split('\t')
        if items[1] == args.p and (args.male == (items[3] == 'male')):
            female_columns[items[0]] = None
#print(female_columns)

#with gzip.open(genome) as fg:
fv = None # variation file handler
skipping = False
with gzip.open(genome) as fg, gzip.open(args.o, 'w') as fo:
    for line in fg:
        if line.startswith('>'):
            if fv is not None: fv.close()
            gpos = 0
            sequence = ''
            chrom = re.split('\\s+', line[1:-1])[0]
            if chrom.startswith('chr') is False:
                chrom = 'chr' + chrom
            if not args.male and chrom.startswith('chrY'):
                skipping = True
                continue
            else:
                skipping = False
            mgen = {}
            fo.write('>{}\n'.format(chrom))
            if chrom not in varfiles:
                sys.stderr.write('missing {}\n'.format(chrom))
                fv = None
            else:
                fv = gzip.open(varfiles[chrom])
                column_determinied = False
                while 1:
                    line = fv.readline()
                    if line == '': break
                    if line.startswith('#CHROM'):
                        items = line.strip().split('\t')
                        for cn, item in enumerate(items):
                            if item in female_columns:
                                female_columns[item] = cn
#                                print(item, cn)
                        column_determined = True
                        break
                if not column_determined:
                    sys.stderr.write('cannot assign columns in {}\n'.format(varfiles[chrom]))
                    fv.close()
                    fv = None
                else:
                    vpos = -1
        else:
            if fv is None or skipping: continue
            # read VCF and cache genotypes
            if vpos < gpos + gmargin:
                mgen_ = {}
                for p_, g_ in mgen.items():
                    if p_ >= gpos: mgen_[p_] = g_
                mgen = mgen_
                while 1:
                    line = fv.readline()
                    if line == '': break
                    items = line.strip().split('\t')
                    pos = int(items[1])
                    rsid = items[2]
                    ref = items[3]
                    alt = items[4]
                    if rsid in preset_snps:
                        b = preset_snps[rsid]
                        if b is None:
                            mgen[pos] = alt
                        else:
                            mgen[pos] = preset_snps[rsid]
                    elif len(ref) == len(alt):
                        genotypes = []
                        for col in female_columns.values():
                            gtype = items[col]
                            if gtype == '0|0':
                                gn = 0
                            elif gtype == '0|1':
                                gn = 1
                            elif gtype == '1|1':
                                gn = 2
                            else:
                                gn = -1
                            genotypes.append(gn)
                        gweight = float(sum(genotypes)) / len(genotypes)
                        replaced = alt if random.random() > gweight else ref
                        if replaced != ref and len(replaced) == 1:
                            mgen[pos] = ref, replaced
                    vpos = pos
                    # if args.verbose and vpos > gpos + vmargin:
                    #     sys.stderr.write('variation file shifted , {} // {}, {}      \n'.format(len(mgen), vpos, gpos  ))
                    #     break
            line = fg.readline().strip()
            for i, base in enumerate(line):
                p_ = gpos + i
                if p_ in mgen:
                    ref, alt = mgen[p_]
                    sys.stderr.write('substituting at {}:{} : {}=>{}        \r'.format(chrom, p_, ref, alt))# base))
                    if args.lowercase:
                        fo.write(alt.lower() if base == ref else complementary[alt].lower())
                    else:
                        fo.write(alt if base == ref else complementary[alt])
                else:
                    fo.write(base)
            fo.write('\n')
            gpos += len(line)
