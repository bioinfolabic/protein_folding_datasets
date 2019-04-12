# -*- coding: utf-8 -*-


def converte(amino):
    if amino == 'G':
        return 'A'
    elif amino == 'A':
        return 'A'
    elif amino == 'L':
        return 'A'
    elif amino == 'M':
        return 'A'
    elif amino == 'F':
        return 'A'
    elif amino == 'W':
        return 'A'
    elif amino == 'K':
        return 'B'
    elif amino == 'Q':
        return 'B'
    elif amino == 'E':
        return 'B'
    elif amino == 'S':
        return 'B'
    elif amino == 'P':
        return 'A'
    elif amino == 'V':
        return 'A'
    elif amino == 'I':
        return 'A'
    elif amino == 'C':
        return 'A'
    elif amino == 'Y':
        return 'B'
    elif amino == 'H':
        return 'B'
    elif amino == 'R':
        return 'B'
    elif amino == 'N':
        return 'B'
    elif amino == 'D':
        return 'B'
    elif amino == 'T':
        return 'B'
    else:
        print 'ERROR 404 amino not found !!!!!!!!!!!!!!!!!!'
        return 0
    
def main(sequencia):
    sequencia_AB = []
    
    sequencia = list(sequencia)

    tamanho = len(sequencia)
    print tamanho
    
    for i in xrange(0,tamanho):
        if sequencia_AB.append(converte(sequencia[i])) == 0 :
            break

    sequencia_AB = ''.join(sequencia_AB)
    print sequencia_AB

    return 'O.K :)'



s = 'GPLGSQRKEKDFQGMLEYHKEDEALLIRNLVTDLKPQMLSGTVPCLPAYILYMCIRHADYTNDDLKVHSLLTSTINGIKKVLKKHNDDFEMTSFWLSNTCRLLHCLKQYSGDEGFMTQNTAKQNEHCLKNFDLTEYRQVLSDLSIQIYQQLIKIAEGVLQPMIVSAMLENESIQGLSGVKPTGYRKRSSSMADGDNSYCLEAIIRQMNAFHTVMCDQGLDPEIILQVFKQLFYMINAVTLNNLLLRKDVCSWSTGMQLRYNISQLEEWLRGRNLHQSGAVQTMEPLIQAAQLLQLKKKTQEDAEAICSLCTSLSTQQIVKILNLYTPLNEFEERVTVAFIRTIQAQLQERNDPQQLLLDAKHMFPVLFPFNPSSLTMDSIHIPACLNLEFLNEV'
print main(s)