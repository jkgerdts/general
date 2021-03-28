
def hist_smoother(csv_path,squishfactor=50):

    import pandas as pd
    from pathlib import Path
    import numpy as np
    
    df = pd.read_csv(csv_path)
    columns = [x for x in df.columns if x not in ['bin','bins']]
    figpath = str(Path(csv_path).parent)+'/'

    def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
            print('error too many dimensions')
            return None
        if x.size < window_len:
            print('error window len too big')
            return None
        if window_len<3:
            return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
            return None
        s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
        #print(len(s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='valid')
        return y

    bins_orig = df['bins'].values
    newbins = [bins_orig[i] for i in range(0,len(bins_orig),squishfactor)]
    df['newbins'] = [next((i for i,v in enumerate(newbins) if v>=x), max(newbins)) for x in bins_orig]
    for col in columns:
        s = sum(df[col].values)
        df[col] = [x/s for x in df[col].values]

    import matplotlib.pyplot as plt

    dfg = df.groupby('newbins').sum()

    for col in columns:
        f,ax = plt.subplots()
        plt.plot(smooth(dfg[col].values,window_len=5))
        plt.title(col)
        plt.savefig(figpath+col+'.pdf')

def get_channels(parent_dir):
    import os
    paths = []
    indexes = []
    for p in os.scandir(parent_dir):
        if p.is_dir():
            paths.append(p.path)
        if "Index.idx.xml" in p.name:
            indexes.append(p.path)
    for pp in paths:
        for p in os.scandir(pp):
            if p.is_dir():
                paths.append(p.path)
            if "Index.idx.xml" in p.name:
                if p.name[0] != '.':
                    indexes.append(p.path)

    import re
    scan = re.compile("Channel: (\d), ChannelName: (DAPI|EGFP|Alexa 647|Alexa 568|Brightfield), Version")

    imnamevars = {'DAPI': 'imname_blue', 'Alexa 647':'imname_farred','Alexa 568':'imname_red','EGFP':'imname_grn','Brightfield':'imname_bf'}

    for ind in indexes:
        print(ind)
        file = open(ind)
        filestr = file.read()
        file.close()
        matches = re.findall(scan,filestr)
        for tup in matches:
            print(imnamevars[tup[1]]+' = "r01c"+zeropad(c)+"f"+zeropad(f)+"p01-ch'+str(tup[0])+'sk1fk1fl1.tiff";')


def makeVirus():
    path = '/Users/jkgerdts/Google Drive/Lab/'
    import subprocess, shutil
    from datetime import datetime
    newfname = datetime.now().strftime("%Y %m_%d")+" virus prep.xlsx"
    vsheet = '293T lentivirus transfection (FuGENE HD, Lipofectamine LTX, GeneJuice).xlsx'
    shutil.copyfile(path+vsheet,path+newfname)
    subprocess.call(('open',path+newfname))

def address():
    print("UCSF Genentech Hall - Room N414")
    print("600 16th st MC2240")
    print("San Francisco CA 94143")

def cloning():
    def openfolder(path):
        from subprocess import call
        call(["open",path])
    openfolder("/Users/jkgerdts/Google Drive/Lab/cloning/")

def pcr():
    print('not yet developed')

def primerTable(df):
    import pandas as pd
    from Bio.SeqUtils.MeltingTemp import Tm_NN as tm
    from Bio.SeqUtils import molecular_weight as mw
    df['mw']=df['seq'].apply(mw)
    df['Tm']=df['seq'].apply(tm)
    return df

def convert(fromunits,tounits,*argv):
    # I'm super proud of this. It converts 
    # handles the following (including prefixed versions)
    # g/l<-->mol/l (g/mol = required argv)
    # g<-->mol (g/mol)
    # X<-->X (note: if scalar entered for tounits the result is a dilution factor)
    # g or mol<-->l (g/mol)(mol/l) OR (g/l) ** detect first argv to decide
    #### would be good to consolidate a bit and remove redundancy; a lot of code is repeated in different if blocks

    def extract(units):
        # this function does the heavy lifting  - parses input strings to separate scalar values from units, keeping track of numerator/denominator and magnitude prefixes
        # returns tuple containing 1) scalar value, 2) measurement type, 3) magintude factor
        import re
        magnitudes={'f':10**-15,'p':10**-12,'n':10**-9,'u':10**-6,'m':10**-3,'c':10**-2,'d':0.1,'k':1000,'':1}
        r=re.compile('([0-9\.e\-]*)([fpnumcdk]?(?!ol)) ?(g|mol|l?)/?([fpnumcdk]?(?!ol)) ?(g|mol|l?)')
        m=r.findall(units)
        scalar = m[0][0]
        if scalar=='':
            scalar=1
        else:
            scalar = float(scalar)
        factor = magnitudes[m[0][1]]/magnitudes[m[0][3]]
        meastype = m[0][2]+'_'+m[0][4]
        #print(scalar,meastype,factor)
        return (scalar,meastype,factor)
    
    def extract_vol(units):
        # this is a much simpler function to extract volume units from a volume only. Should return a string like "ml"
        import re
        try:
            r=re.compile('[fpnumcdk]?l')
            m=r.findall(units)
            return m[0]
        except:
            print(units,m)
        return None
        

    (f_s,f_m,f_f)=extract(fromunits)
    (t_s,t_m,t_f)=extract(tounits)
    
    conversion = f_m+'-->'+t_m
    factor = 0
    printstring = 'this was not set'
    
    if f_m==t_m:
        factor = (f_s*f_f/(t_s*t_f))
        if t_s==1:
            printstring = fromunits+' is @ '+tounits
        else:
            printstring = 'dilution factor is @'
    elif conversion in ['g_l-->mol_l','g_-->mol_']:
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='g_mol':
            factor = (f_s*f_f)/(a_s*a_f)/(t_s*t_f)
            printstring = printstring = fromunits+' is @ '+tounits
        else:
            print('error: expected single argv entry Xg/mol but received '+a_m)
            
    elif conversion in ['g_l-->l_']:
        (a_s,a_m,a_f)=extract(argv[0])
        volunits = extract_vol(tounits)
        if a_m=='g_l':
            factor = (f_s*f_f)/(a_s*a_f)
            volume = t_s/factor
            volume_dil = t_s-volume
            printstring = 'dilution factor is @; mix '+"{:.3f}".format(volume)+' '+volunits+' into '+"{:.3f}".format(volume_dil)+' '+volunits+' solvent'
        elif a_m=='g_mol':

            ############################
            # next argv must be concentration in mol_l
            (a2_s,a2_m,a2_f)=extract(argv[1])
            
            if a2_m=='mol_l':
            
                factor = (f_s*f_f)/(a_s*a_f)/(a2_s*a2_f)
                volume = t_s/factor
                volume_dil = t_s-volume
                #  should be g/l * mol/g * l/mol
                factor_vol = factor*(t_s*t_f) # this is a quantity in L and I don't know what it means
                print('factor_vol = '+"{:.3f}".format(factor_vol)+' l (?)')
                printstring = 'to dilute a solution of '+fromunits+' to '+tounits+' of '+argv[1]+' given molar mass '+argv[0]+', dilution factor is @; mix '+"{:.3f}".format(volume)+' '+volunits+' into '+"{:.3f}".format(volume_dil)+' '+volunits+' solvent'
            
            else:
                print('error: expected 4th argument to be mol_l; not sure what to do with current input')
        
            ###########################
            
            
        else:
            print('error: not sure what to do with these units')

    elif conversion in ['mol_l-->l_']:
        (a_s,a_m,a_f)=extract(argv[0])
        volunits = extract_vol(tounits)
        if a_m=='mol_l':
            factor = (f_s*f_f)/(a_s*a_f)
            volume = t_s/factor
            volume_dil = t_s-volume
            printstring = 'dilution factor is @; mix '+"{:.3f}".format(volume)+' '+volunits+' into '+"{:.3f}".format(volume_dil)+' '+volunits+' solvent'
        elif a_m=='g_mol':
            
            ############################
            # next argv must be concentration in g_l
            (a2_s,a2_m,a2_f)=extract(argv[1])
            
            if a2_m=='g_l':
            
                factor = (f_s*f_f)*(a_s*a_f)/(a2_s*a2_f)
                volume = t_s/factor
                volume_dil = t_s-volume
                #  should be mol/l * g/mol * l/g
                factor_vol = factor*(t_s*t_f) # this is a quantity in L and I don't know what it means
                print('factor_vol = '+"{:.3f}".format(factor_vol)+' l (?)')
                printstring = 'to dilute a solution of '+fromunits+' to '+tounits+' of '+argv[1]+' given molar mass '+argv[0]+', dilution factor is @; mix '+"{:.3f}".format(volume)+' '+volunits+' into '+"{:.3f}".format(volume_dil)+' '+volunits+' solvent'
            
            else:
                print('error: expected 4th argument to be g_l; not sure what to do with current input')
        
            ###########################
        else:
            print('error: not sure what to do with these units')
            
    elif conversion in ['mol_l-->g_l','mol_-->g_']:
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='g_mol':
            factor = (f_s*f_f)*(a_s*a_f)/(t_s*t_f)
            printstring = printstring = fromunits+' is @ '+tounits
        else:
            print('error: expected single argv entry Xg/mol but received '+a_m)
    elif conversion=='g_-->l_':
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='g_l':
            factor = (f_s*f_f)/(a_s*a_f)/(t_s*t_f)
            printstring = 'mix '+fromunits+' into @ '+tounits+' to make '+argv[0]
        else:
            (a2_s,a2_m,a2_f)=extract(argv[1])
            if a_m+'__'+a2_m=='g_mol__mol_l':
                factor = (f_s*f_f)/(a_s*a_f)/(a2_s*a2_f)/(t_s*t_f)
                printstring = 'mix '+fromunits+' into @ '+tounits+' to make '+argv[1]
            else:
                print('error, anticipated g/mol then mol/l argvs')
                
    elif conversion=='mol_-->l_':
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='mol_l':
            factor = (f_s*f_f)/(a_s*a_f)/(t_s*t_f)
            printstring = 'mix '+fromunits+' into @ '+tounits+' to make '+argv[0]
        else:
            (a2_s,a2_m,a2_f)=extract(argv[1])
            if a_m+'__'+a2_m=='g_mol__g_l':
                factor = (f_s*f_f)*(a_s*a_f)/(a2_s*a2_f)/(t_s*t_f)
                printstring = 'mix '+fromunits+' into @ '+tounits+' to make '+argv[1]
            else:
                print('error, anticipated g/mol then g/l argvs')

    elif conversion=='l_-->g_':
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='g_l':
            factor = (f_s*f_f)*(a_s*a_f)/(t_s*t_f)
            printstring = 'mix @ '+tounits+' into '+fromunits+' to make '+argv[0]
        else:
            (a2_s,a2_m,a2_f)=extract(argv[1])
            if a_m+'__'+a2_m=='g_mol__mol_l':
                factor = (f_s*f_f)*(a2_s*a2_f)*(a_s*a_f)/(t_s*t_f)
                printstring = 'mix @ '+tounits+' into '+fromunits+' to make '+argv[1]
            else:
                print('error, anticipated g/mol then mol/l argvs')
                
    elif conversion=='l_-->mol_':
        (a_s,a_m,a_f)=extract(argv[0])
        if a_m=='mol_l':
            factor = (f_s*f_f)*(a_s*a_f)/(t_s*t_f)
            printstring = 'mix @ '+tounits+' into '+fromunits+' to make '+argv[0]
        else:
            (a2_s,a2_m,a2_f)=extract(argv[1])
            if a_m+'__'+a2_m=='g_mol__g_l':
                factor = (f_s*f_f)*(a2_s*a2_f)/(a_s*a_f)/(t_s*t_f)
                printstring = 'mix @ '+tounits+' into '+fromunits+' to make '+argv[0]
            else:
                print('error, anticipated g/mol then g/l argvs')
                printstring = 'mix @ '+tounits+' into '+fromunits+' to make '+argv[1]

    else:
        print('not yet implemented')
    if any([factor > 1000,factor < 0.001]):
        factor_str = format(factor,'.2e')
    elif factor > 100:
        factor_str = format(factor,'.0f')
    elif factor > 10:
        factor_str = format(factor,'.1f')
    elif factor  > 1:
        factor_str = format(factor,'.2f')
    elif factor  > 0.1:
        factor_str = format(factor,'.3f')
    elif factor  > 0.01:
        factor_str = format(factor,'.4f')
    else:
        factor_str = format(factor,'.5f')

    


    print(printstring.replace('@',factor_str))
    return factor

    

def getJWPlasmids():
    import os
    root = '/Users/jkgerdts/Desktop/Plasmids'
    ff=[]
    for root, dirs, files in os.walk('/Users/jkgerdts/Desktop/Plasmids'):
        for name in files:
            file = os.path.join(root, name)
            if '.dna' in file:
                ff.append(file)
    return ff

def JWBlastn(st):
    files = getJWPlasmids()
    from snapgene_reader import snapgene_file_to_dict as SG
    for ff in files:
        s=SG(ff)
        if st.lower() in s['seq'].lower():
            print(ff)
            
        
def getLenths(fnames):
    import numpy as np
    from snapgene_reader import snapgene_file_to_dict as SG
    arr = np.zeros(len(fnames))
    for i,f in enumerate(fnames):
        s=SG(f)
        arr[i]=len(s['seq'])
    for percentile in [10,25,50,75,90,95]:
        print(percentile,np.percentile(arr,percentile))
    return arr




def syncodons(s,npadding,enz):


    default_enzymes = ['AbsI', 'Acc65I', 'AjuI', 'AsiSI', 'BclI', 'BlpI', 'BmgBI', 'BsiWI', 'BstBI', 'BstXI', 'BstZ17I', 'Bsu36I', 'FseI', 'FspAI', 'KpnI', 'MauBI', 'MluI','MreI', 'MscI', 'NdeI', 'PacI', 'PasI', 'PmeI', 'PspXI', 'PsrI', 'RsrII', 'SgrDI', 'SnaBI', 'SpeI', 'SrfI', 'SwaI']
    if len(enz)<1:
        enz=default_enzymes


    #seq2a = 'AGGAAGCGACGGAGTGGGTCAGGCGAGGGGAGAGGCTCACTCCTCACGTGTGGAGACGTTGAAGAAAACCCCGGACCA'
    #test = 'AGGAAGCGA'
    #s='nnnnnnnnn'+seq2a
    s = 'n'*npadding+s
    codons = [s[i:i+3] for i in range(0,len(s),3)]

    syncodons = []

    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    rgencode = {}
    for v in gencode.values():
        rgencode[v]=[]
    for k,v in gencode.items():
        rgencode[v].append(k)

    import re
    for codon in codons:
        if 'n' in codon:
            # these will be handled without translation table
            pcodons = [list(gencode.keys())[i] for i,v in enumerate([re.fullmatch(codon.replace('n','[GATC]'),x) is not None for x in gencode.keys()]) if v]
            additions = []
            for c in pcodons:
                synonyms = rgencode[gencode[c]]
                for syn in synonyms:
                    if syn not in pcodons:
                        additions.append(syn)
            for syn in additions:
                pcodons.append(syn)
        else:
            # need to check translation table
            pcodons = rgencode[gencode[codon]]
        syncodons.append(pcodons)
    
    old = ['']
    new = []
    for i, sc in enumerate(syncodons):
        for scc in sc:
            for oldstr in old:
                new.append(oldstr+scc)
        old=new.copy()
        new = []
    #return syncodons
    synstrings = old.copy()

    def getrestrictionmatches(seq,starting,enz):
        from Bio import Restriction
        rb = Restriction.RestrictionBatch(enz)
        from Bio import Seq
        from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
        amb = IUPACAmbiguousDNA()
        s = Seq.Seq(seq,amb)
        dic = rb.search(s)
        hits = []
        sites = []
        for k,v in dic.items():
            if len(v)>0:
                if max(v)>=starting:
                    hits.append(k)
                    sites.append(max(v))
        return (hits,sites)

    hitsdic = {}

    for syn in synstrings:
        (hits,sites) = getrestrictionmatches(syn,npadding,enz)
        for i, h in enumerate(hits):
            newsite = True
            if h in hitsdic.keys():
                if sites[i] in hitsdic[h]:
                    newsite = False
                else:
                    hitsdic[h].append(sites[i])
            else:
                hitsdic[h]=[sites[i]]
            if newsite:
                print(h,sites[i],syn)




def randomizecodons(s,exclude_rare=True):
    rarecodons = ['CGA','CGG','TCG','CGC','CCG','GCG','ACG','CGT']
    import numpy as np
    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    rgencode = {}
    for v in gencode.values():
        rgencode[v]=[]
    for k,v in gencode.items():
        rgencode[v].append(k)
    codons = [s[i:i+3] for i in range(0,len(s),3)]
    syncodons = []
    for c in codons:
        pcodons = [cod for cod in rgencode[gencode[c]] if cod not in rarecodons]
        syncodons.append(np.random.choice(pcodons))
    return ''.join(syncodons)


def quintara():
    print('SFFV Fwd    QB3969	AATGACCCTGCGCCTTATTT')
    print('WPRE-Rev	QB0112	CATAGCGTAAAAGGAGCAACA')
    print('IRES-For	QB0037	TGGCTCTCCTCAAGCGTATT')
    print('IRES-Rev	QB0038	CCTCACATTGCCAAAAGACG')
    print('Rev.mCherry(ON)	QB2307	CCTCGATCTCGAACTCGTGGC')


def read_flojo(df):
    #Specimen_011_B8_B08_031.fcs
    #Specimen_011_B8_B08_031.fcs/live cells/K562/negative/activated
    import re
    parser = re.compile('([A-Za-z0-9_\s]+)\.fcs\/([A-Za-z0-9_\s]+)?\/([A-Za-z0-9_\s]+)?\/([A-Za-z0-9_\s]+)?\/([A-Za-z0-9_\s]+)?\/([A-Za-z0-9_\s]+)?')
    

def rch2num(r):
    return ord(r.upper())-64

def rnum2ch(r):
    return chr(r+64)

def well2rc(well,plate=96):
    rowlen = {6:3,12:4,24:6,96:12,384:24}
    r = int((well-1)/rowlen[plate])+1
    c = well - (r-1)*rowlen[plate]
    return (r,c)

def well2rcstr(well,plate=96,rcfmt = 'R - 0C'):
    (r,c) = well2rc(well,plate)
    rcfmt = rcfmt.replace('0C','{:0>2d}'.format(c))
    rcfmt = rcfmt.replace('C',str(c))
    rcfmt = rcfmt.replace('R',rnum2ch(r))
    return rcfmt  

def rcstr2well(rcstr,plate=96):
    import re
    parser = re.compile('([A-Pa-p])\s?-?_?\s?([0-9]{1,2})')
    parsed = parser.findall(rcstr)

    rowlen = {6:3,12:4,24:6,96:12,384:24}
    try:
        r = rch2num(parsed[0][0].upper())
        c = int(parsed[0][1])
        well = (r-1)*rowlen[plate]+c
        #print(r,c,well)
    except:
        print('parsing error -jkg')
        return 0
    return well


def reviewboards():
    
    input('copy the whole table and press enter')
    
    import pandas as pd
    import random

    df = pd.read_clipboard()

    def shuffle(dff):
        l = len(dff)
        arr = [x for x in range(l) if dff['W'].loc[x]*0.1>random.random()]
        random.shuffle(arr)
        it = iter(arr)
        print('shuffling; returning: ',len(arr))
        return it

    it = shuffle(df)
    wcol = df.columns.get_loc('W')

    resp = ''
    counter = 0
    while resp not in ['X']:
        if counter > 10:
            it = shuffle(df)
            counter = 0
        counter = counter+1
        i = next(it)
        print('QUESTION',df.iloc[i]['S'])
        resp = input('')
        print('ANSWER',df.iloc[i]['R'])
        weight = input('')
        try:
            df.iloc[i,wcol]=int(weight)
        except:
            print('error assigning weight from input',weight)
        print('----------------------')
    
    df.to_clipboard()
    print('dataframe copied to clipboard and returned from function call')
    return df
        
