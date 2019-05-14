# Pavithra Nagarajan
# May 2019
# Purpose of this Project: Generate optimal RNA construct with randomly generated buffer and hairpin regions that do NOT disrupt native RNA folding.
# In order to validate no structural variance, heavily relied on Matthews Lab Web Server for RNA secondary structure visualization.
# Designed to be run BY rna_construct.m !!!

# In essence, this script generates a text file with RNA construct candidates.
# These candidates are most likely non-disruptive, but the MATLAB file actually runs the RNAStructure package to ensure that the RNA constructs are valid!

import itertools
import re
import random


class CommandLine():
    '''
    Handles user input arguments.
    '''

    def __init__(self, inOpts=None):

        import argparse
        self.parser = argparse.ArgumentParser(description='user specifications for rna_construct.py',
                                              epilog='This tool will help you generate a (hopefully optimal) RNA construct with buffer and hairpin regions.',
                                              add_help=True,  # default is True
                                              prefix_chars='-',
                                              usage='%(prog)s [options] -option1[default] <input >output'
                                              )
        self.parser.add_argument('outputfile', action='store', help='output file name')
        self.parser.add_argument('rna', action='store', help='rna sequence')
        self.parser.add_argument('-buffer1', '--buffer1', action='store', nargs='?', const=True, default=False, help='buffer1 dot bracket check')
        self.parser.add_argument('-buffer2', '--buffer2', action='store', nargs='?', const=True, default=False, help='buffer2 dot bracket check')
        self.parser.add_argument('-buffer3', '--buffer3', action='store', nargs='?', const=True, default=False, help='buffer3 dot bracket check')  # allows multiple list options
        self.parser.add_argument('-hairpin', '--hairpin', action='store', nargs='?', const=True, default=False, help='hairpin dot bracket check')
        self.parser.add_argument('-hairpincheck', '--hairpincheck', action='store', nargs='?', const=True, default=False, help='hairpin partial complement check')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def generate_complement1(sequence):
    '''Simply translates passed in RNA sequence into its canonical complement. Does not reverse direction. '''
    sequence = sequence.upper()
    invalid = re.search("[^ACUG]", sequence)
    if invalid:
        print('Error: Sequence passed in to generate_complement() is not a valid RNA sequence.\n')
        print("{0} detected.\n".format(invalid.group()))
        sys.exit()

    sequence = ''.join(sequence).translate(str.maketrans('ACGU', 'UGCA'))
    return(sequence)


def find_substrings(sequence):
    '''Finds overlapping substrings of size 3."""
            Example:
            Input: 'ACTGTGA'
            Output: ['ACT', 'CTG', 'TGT', 'TGA']
    '''
    pattern = '(?=(.{3}))'  # Lookahead assertion "?="" enables overlapping search
    result = re.findall(pattern, sequence)
    return result


def create_hairpin(length=5, loop_length=5):
    '''Creates a hairpin sequence, of default size 5 for the helix, and 5 for the loop.
    '''

    # Create hairpin helix
    # 2 are fixed to be GC
    one_side = ''.join([random.choice('ACUG') for i in range(length - 2)]) + 'CG'
    # Checks to see if any consecutive runs found. If so, alternate to make more stable.
    c = 0
    check = False
    execution_count = 0
    # Finds distinct groups
    for group, values in itertools.groupby(one_side):
        v = list(values)
        if check == True:
            break
        for i in v:
            if len(v) > 1:
                bases = v
                c += 1
                check = True
                break

    # if consecutive identical nucleotides have been found in hairpin helix...
    while(c != 0):

        # To prevent infinite loop...
        if execution_count > 2:
            one_side = ''.join([random.choice('ACUG') for i in range(length - 2)]) + 'CG'
            execution_count = 0
            bases = []

        if len(bases) > 1:
            # i.e. 'AAA' or 'UU'
            pattern = ''.join(bases)

            # If consecutive fragment has even length
            if len(pattern) % 2 == 0:
                # Search for the position of that
                find = re.search(pattern, one_side)
                # Store the starting index of the fragment
                i = find.span()[0]
                count = 0

                while count < (len(pattern) / 2):
                    one_side = list(one_side)
                    # Switch the nucleotide
                    one_side[i] = generate_complement1(one_side[i])
                    count += 1
                    i += 2
                    one_side = ''.join(one_side)

            # If consecutive fragment has odd length
            if len(bases) % 2 != 0:
                # Search for the position of it
                find = re.search(pattern, one_side)
                i = find.span()[0]
                count = 0
        # if it's 3... need to loop twice 3 divided by 2 is 1
        # if it's 5...need to loop 3 times. 5 divided by 2 is 1
                while i < (len(bases) // 2):
                    one_side = list(one_side)
                    count += 1
                    one_side[i + 1] = generate_complement1(one_side[i + 1])
                    i += 2
                    one_side = ''.join(one_side)

        # Reset
        c = 0
        check = False
        # Finds distinct groups
        for group, values in itertools.groupby(one_side):
            v = list(values)
            if check == True:
                break
            for i in v:
                if len(v) > 1:
                    bases = v
                    c += 1
                    check = True
                    break
        execution_count += 1

    # Creating Hairpin helix's canonical reverse complement
    complement = ''.join((one_side).translate(str.maketrans('ACGU', 'UGCA'))[::-1])

    # Create the Loop, of size 5
    loop = ''.join([random.choice('ACUG') for i in range(loop_length)])
    # Ensuring that loop will stay as a loop, and not base pair with nucleotides across it.
    # Accounting for GU wobble pairing as well, which is why 2 dictionaries needed.
    pairing1 = {'A': 'U', 'G': 'U', 'C': 'G', 'U': 'A'}
    pairing2 = {'A': 'U', 'G': 'C', 'C': 'G', 'U': 'G'}
    i = 0
    while i < (loop_length / 2) and (pairing1[(loop[i])] == loop[loop_length - (i + 1)] or pairing2[loop[i]] == loop[loop_length - (i + 1)]):
        loop = ''.join([random.choice('ACUG') for i in range(length)])
        i += 1

    # Now, time to ensure loop is not complementary to rest of hairpin. It's not enough to check that it's not complementary with itself.
    loop_substrings = find_substrings(loop)
    check = True
    i = 0
    while i < len(loop_substrings):
        for x in generate_complement3(one_side + complement):
            check = re.findall('(' + loop_substrings[i] + ')', x[::-1])  # check in reverse
            if check:
                i = 0
                loop = ''.join([random.choice('ACUG') for i in range(length)])  # regenerate loop
                loop_substrings = find_substrings(loop)  # regenerate substrings
        i += 1
    # Combine Hairpin Helix and Loop together, to form final hairpin structure.
    hairpin = one_side + loop + complement
    return (hairpin)


def generate_complement2(seq):
    '''Generates all possible combinations of compliments (both canonical, and non-canonical UG pairing) for a given 3-nucleotide RNA fragment.
       Assumption: ONLY SEQUENCE FRAGMENTS OF LENGTH 3.
        '''
    G_positions = [index.start() for index in re.finditer('G', seq)]
    U_positions = [index.start() for index in re.finditer('U', seq)]
    positions = G_positions + U_positions

    compliments = {'G': 'U', 'U': 'G'}
    compliments_reg = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

    compliment = set()
    seq = list(seq)
    original = ''.join(list(seq))

    combos = list(map(list, list(itertools.combinations(positions, 1)) + list(itertools.combinations(positions, 2))))
    for mini_list in combos:
        for base in mini_list:
            seq[base] = compliments[seq[base]]

        # Set difference
        if len(mini_list) > 1:
            s2 = list({0, 1} - set(mini_list))
            for base in s2:
                seq[base] = compliments_reg[seq[base]]

        # Add the possible compliment to set of compliments
        compliment.add(''.join(seq))

        # Reset
        seq = list(original)

    compliment = list(compliment)

    # Make sure to include basic, canonical RNA compliment
    compliment.append(''.join(seq).translate(str.maketrans('ACGU', 'UGCA')))
    return list(compliment)


def generate_complement3(seq):
    '''Generates all possible combinations of compliments (both canonical, and non-canonical UG pairing) for a given 3-nucleotide RNA fragment.
       Assumption: ONLY SEQUENCE FRAGMENTS OF LENGTH 3.
        '''
    G_positions = [index.start() for index in re.finditer('G', seq)]
    U_positions = [index.start() for index in re.finditer('U', seq)]
    positions = G_positions + U_positions

    compliments = {'G': 'U', 'U': 'G'}
    compliments_reg = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

    compliment = set()
    seq = list(seq)
    original = ''.join(list(seq))

    combos = list(map(list, list(itertools.combinations(positions, 1)) + list(itertools.combinations(positions, 2)) + list(itertools.combinations(positions, 3))))
    for mini_list in combos:
        for base in mini_list:
            seq[base] = compliments[seq[base]]

        # Set difference
        if len(mini_list) > 1:
            s2 = list({0, 1} - set(mini_list))
            for base in s2:
                seq[base] = compliments_reg[seq[base]]

        # Add the possible compliment to set of compliments
        compliment.add(''.join(seq))

        # Reset
        seq = list(original)

    compliment = list(compliment)

    # Make sure to include basic, canonical RNA compliment
    compliment.append(''.join(seq).translate(str.maketrans('ACGU', 'UGCA')))
    return list(compliment)


def check_buffer(buffer_of_interest, reference_buffer):
    complements1 = []
    comps1 = []
    for i in range(len(reference_buffer)):
        x = generate_complement2(reference_buffer[i:i + 2])
        comps1.append(x)
        for c in range(len(x)):
            c = frozenset([x[c]])
            complements1.append(c)
    complements1 = set(complements1)
    complements1 = list(map(''.join, complements1))
    check = True
    for x in complements1:
        if x in buffer_of_interest:
            check = False
    return check


def check_hairpin(buffer1, buffer2, buffer3, hairpin_check=10):
    complements1 = []
    comps1 = []
    for i in range(len(buffer1)):
        x = generate_complement3(buffer1[i:i + 3])
        comps1.append(x)
        for c in range(len(x)):
            c = frozenset([x[c]])
            complements1.append(c)
    complements1 = set(complements1)
    complements1 = list(map(''.join, complements1))
    complements2 = []
    comps2 = []
    for i in range(len(buffer2)):
        x = generate_complement3(buffer2[i:i + 3])
        comps2.append(x)
        for c in range(len(x)):
            c = frozenset([x[c]])
            complements2.append(c)
    complements2 = set(complements2)
    complements2 = list(map(''.join, complements2))

    check = False
    tracker = 0
    while(check == False) and tracker < 100:
        tracker += 1
        hairpin = create_hairpin()
        check = True
        for x in complements2:
            if x in hairpin:
                check = False
        if check == True:
            count = 0
            for y in complements1:
                if y in hairpin:
                    count += 1
            if count < hairpin_check:
                return hairpin
    # print('No suitable hairpin found with the given input for hairpin_check. Please increase value.\n')
    return ''


def main(inCL=None):

    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    file = open(myCommandLine.args.outputfile, "w")
    rna = myCommandLine.args.rna

    # All 3 length complements of RNA, also including UG base pairing
    complements = []
    comps = []
    for i in range(len(rna)):
        x = generate_complement3(rna[i:i + 3])
        comps.append(x)
        for c in range(len(x)):
            c = frozenset([x[c]])
            complements.append(c)
    complements = set(complements)
    complements = list(map(''.join, complements))

    # All possible 10-length buffer candidates that are not-complementary to RNA
    c = list(map(''.join, list(itertools.product(['A', 'C', 'U', 'G'], repeat=10))))
    [frozenset([x]) for x in c]
    c = list(set(c))
    possibles = []
    counter = 0
    i = 0
    count = 0
    for i in c:
        count = 0
        for x in complements:
            if x in i:
                count += 1
        if count < 2:
            possibles.append(i)

    # All possible choices of 3 buffers, with identical buffers allowed.
    buffers = list(itertools.product(possibles, repeat=3))
    # Prints length, helps to get an idea of whether this RNA is going to be tough to generate regions that don't disrupt its natural folding, or not.
    file.write("Number of choices for 3 buffers, including duplicates: {0}\n".format(len(buffers)))

    for a in range(len(buffers)):
        (buffer1, buffer2, buffer3) = buffers[a]

        if check_buffer(buffer2, buffer1) and check_buffer(buffer3, buffer1) and check_buffer(buffer3, buffer2):
            hairpin = check_hairpin(buffer1, buffer2, buffer3, int(myCommandLine.args.hairpincheck))
            if hairpin != '':
                file.write(buffer1 + hairpin + buffer2 + rna + buffer3)
                file.write('\n')
    print("{0}{1}{2}{3}".format(int(myCommandLine.args.buffer1), int(myCommandLine.args.buffer2), int(myCommandLine.args.buffer3), int(myCommandLine.args.hairpin)))


main()


"""
RNA
rna_of_interest = 'CCCGCCUGGAGGCCGCGGUCGGCCCGGGGCUUCUCCGGAGGCACCCACUGCCACCGCGAAGAGUUGGGCUCUGUCAGCCGCGGG';
rna_of_interest='GGCCAAAGGCGUCGAGUAGACGCCAACAACGGGUAUGCACAGAUGUGGAAACAGGAACUGAUGUGUCCAUUACACCACUAGGACAGAGGCCAGAACAAUGAAGAAACCAAAUACUUGGAAGAGGGUAGAGAUAAUGAAUGGAGUCCAAGAGCCCUGAUUGUGCCAUAAAUGUCCAGAUAAUUCCAUACCUGAGGAUUAUGUGGUUUGUAAACUUGGCACUUAGAAGAACCAAUAAAAUCAUGUUAUAGUUUCAAAAACCAAACCGUCAGCGAGUAGCUGACAAAAAGAAACAACAACAACAAC'
rna_of_interest='GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCA'
"""
