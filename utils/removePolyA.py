import sys
import mappy

def removePolyA(seq):
    reverse = seq[::-1]
    Astate, Astretch, Vstretch, trimPos = False, 0, 0, 0
    for pos in range(0, len(reverse), 1):
        base = reverse[pos]
        if not Astate:
            if base == 'A':
                Astretch += 1
                if Astretch == 6:
                    Astate = True
                    lastA = pos
            else:
                Astretch=0
        if Astate:
            if base != 'A':
                Vstretch += 1
                Astretch = 0
            else:
                Astretch += 1
                if Astretch >= 3:
                    Vstretch = 0
                    lastA = pos
            if Vstretch >= 3:
                trimPos = lastA
                break
    reverseTrim = reverse[trimPos:]
    seqTrim = reverseTrim[::-1]
    return seqTrim, Astate


def main():
    trim=False
    if len(sys.argv)>3:
        trim=sys.argv[3].split(',')
        fivePrime=int(trim[0])
        threePrime=-int(trim[1])
    out=open(sys.argv[2],'w')
    for name,seq,qual in mappy.fastx_read(sys.argv[1]):
        if trim:
            seq=seq[fivePrime:threePrime]
        seqTrim,Astate=removePolyA(seq)
        out.write('>'+name+'\n'+seqTrim+'\n')
    out.close()

main()
