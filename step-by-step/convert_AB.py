# -*- coding: utf-8 -*-


def convert(amino):
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
    
def main(sequence, name):
    name = name + '.txt'
    sequence_AB = []
    
    sequence = list(sequence)

    tamanho = len(sequence)
    print tamanho
    
    for i in xrange(0,tamanho):
        if sequence_AB.append(convert(sequence[i])) == 0 :
            break

    sequence_AB = ''.join(sequence_AB)
    print sequence_AB

    arquivo = open(name, 'w')
    arquivo.write(sequence_AB)
    arquivo.close()

    return 'O.K :)'



s = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
name = '2GB1'
print main(s, name)