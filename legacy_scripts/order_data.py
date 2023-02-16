""" After extractor.py is used to obtain a number of data files, this script compiles and unscrambles that data. """

import getopt
import numpy as np
import sys

def usage():
    print('usage')

def main():
    # Get options and arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'l:r:', ['help'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit()

    for o, a in opts:
        if o in ('--help'):
            usage()
            sys.exit()
        elif o in ('-l'):
            remote = ''
            save_name = a
        elif o in ('-r'):
            remote = '-remote'
            save_name = a
        else:
            assert False, "unhandled option"

    try:
        save_name
    except:
        print("Must choose local (-l) or remote (-r).")
        sys.exit()

    save_path = f'saves{remote}/{save_name}'

    # Gather the q_data files into a list
    q_data = []
    for i in range(56):
        with open(f'{save_path}/chk/data_{i}.txt') as file:
            q_data.append(file.read())

    # Write the q_data to one q_data.txt
    with open(f'{save_path}/chk/q_data.txt','w') as file:
        for i in range(56):
            file.write(q_data[i])

    # Get unscramble key as ordered list of positions
    with open(f'{save_path}/chk/scramble-key.txt','r') as file:
        key = file.read().splitlines()
    key = [int(num) for num in key]

    # Dimension of data and length of data
    d = 5
    n = len(key)

    # Open q_data.txt and reshape to numpy array
    with open(f'{save_path}/chk/q_data.txt','r') as file:
        data_unordered = np.fromstring(file.read(),sep=' ').reshape((-1,d))

    # Create blank numpy array for the ordered data
    data_ordered = np.zeros((n,d))

    # Unscramble the unordered data and write into ordered data
    for k in range(n):
        data_ordered[key[k]] += data_unordered[k]

    # Write the ordered data to text file
    with open(f'{save_path}/chk/q_data_unsrambled.txt','w') as file:
        for vec in data_ordered:
            vec_str = [str(coord) for coord in vec]
            file.write(' '.join(vec_str) + '\n')
    
    print('Done writing ordered data.')

if __name__ == '__main__':
    main()