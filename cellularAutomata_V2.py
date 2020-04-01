import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import colors
import re								# Regular expression - used for parsing

# Code is referenced from python sympy library. Reference for further explination


def rref(B, tol=1e-8, debug=False):
    A = B.copy()
    rows, cols = A.shape
    r = 0
    pivots_pos = []
    row_exchanges = np.arange(rows)
    for c in range(cols):
        if debug:
            print("Now at row", r, "and col", c, "with matrix:")
            print(A)

        # Find the pivot row:
        pivot = np.argmax(np.abs(A[r:rows, c])) + r
        m = np.abs(A[pivot, c])
        if debug:
            print("Found pivot", m, "in row", pivot)
        if m <= tol:
            # Skip column c, making sure the approximately zero terms are
            # actually zero.
            A[r:rows, c] = np.zeros(rows-r)
            if debug:
                print("All elements at and below (", r,
                      ",", c, ") are zero.. moving on..")
        else:
            # keep track of bound variables
            pivots_pos.append((r, c))

            if pivot != r:
                # Swap current row and pivot row
                A[[pivot, r], c:cols] = A[[r, pivot], c:cols]
                row_exchanges[[pivot, r]] = row_exchanges[[r, pivot]]

                if debug:
                    print("Swap row", r, "with row", pivot, "Now:")
                    print(A)

            # Normalize pivot row
            A[r, c:cols] = A[r, c:cols] / A[r, c]

            # Eliminate the current column
            v = A[r, c:cols]
            # Above (before row r):
            if r > 0:
                ridx_above = np.arange(r)
                A[ridx_above, c:cols] = A[ridx_above, c:cols] - \
                    np.outer(v, A[ridx_above, c]).T
                if debug:
                    print("Elimination above performed:")
                    print(A)
            # Below (after row r):
            if r < rows-1:
                ridx_below = np.arange(r+1, rows)
                A[ridx_below, c:cols] = A[ridx_below, c:cols] - \
                    np.outer(v, A[ridx_below, c]).T
                if debug:
                    print("Elimination below performed:")
                    print(A)
            r += 1
        # Check if done
        if r == rows:
            break
    return (A, pivots_pos, row_exchanges)


def cellular_automata(num_elements, num_alphabet, cellular_automata, num_steps, evolution_matrix):
    """Takes a state and evolves it over n steps.

    The final output will look like:
            0001
            1001
            0101
            1111
    """

    ca_next = cellular_automata[:]  # List representation of next state

    print(*ca_next)  # Prints the starting state

    step = 1  # Number of current state

    M = []
    M.append(ca_next)

    while (step < num_steps):
        ca_next = []  # Reset list for every new step

        # Add elements to next state according to update rule
        for i in range(0, num_elements):
            if i > 0:
                ca_next.append(
                    (cellular_automata[i-1]+cellular_automata[i]) % num_alphabet)

            elif(i == 0):
                ca_next.append(
                    (cellular_automata[num_elements - 1]+cellular_automata[i]) % num_alphabet)

        M.append(ca_next)

        cellular_automata = ca_next[:]  # Update cell list

        step += 1  # Step increment

    for i in range(0, num_steps):
        print(M[i], end = " ")
        print()

    # plt.matshow(M)
    # plt.show()

    #print(M)

    M_rref = np.asarray(M, dtype=np.int32)
    print(rref(M_rref))
    # rref(M_rref, tol=1e-8, debug=True)


def check_state(num_elements):
    """
    This process takes a string of integers as input and verifies that it is a valid starting state.
    """

    start_state = []  # Starting state will be appended one element at a time.

    num_digits = 0
    test_state = input('Enter starting state: ')
    for i in test_state:
        if (i.isdigit()):
            num_digits = num_digits + 1
        else:
            print('Incorrect character')

    while (num_digits != num_elements):
        print('You entered: ', num_digits,
                ' element(s)\nThis automaton needs: ', num_elements, ' element(s)\n')
        num_digits = 0
        test_state = input('Enter starting state: ')
        for i in test_state:
            if (i.isdigit()):
                num_digits = num_digits + 1
            else:
                print('Incorrect character')

    for i in test_state:
        start_state.append(int(i))
        
    return start_state


def evolve_matrix(num_elements, update_rule, debug=False):
    """
    This function translates an identity_matrix into an evolution_matrix, given an update rule.
        num_elements        -- Number of cells in a row.
        update_rule         -- Rule that defines the evolution_matrix.
        debug               -- Flag that offers insight to program if set.
        identity_matrix     -- n by n matrix filled with 0s aside from 1s in the top-left diagonal.
        evolution_matrix    -- n by n matrix that when multiplied to a state, gives the next state.
        row                 -- Temporary variable for each row in the matrix.
    """

    identity_matrix = np.identity(num_elements, int)
    print(identity_matrix)

    evolution_matrix = []
    row = []
    new_element = 0

    for i in range(0, num_elements):		# Every row in matrix
        row = []
        for j in range(0, num_elements):    # Every element in row
            new_element = 0
            
            for k in update_rule:           # Every element in update rule
                if j + k >= num_elements:
                    l = j + k - num_elements
                else:
                    l = j + k
                    
                if debug == True:
                    print('l = ', l)
                    
                new_element = new_element + identity_matrix[i][l]
                    
                if debug == True:           # Debug -- print the location and value of element
                    print('Element [', i, ', ', j, '] will be ', new_element)
            
            row.append(new_element % num_alphabet)
            
            if debug == True:
                print(new_element, ' % ', num_alphabet, ' = ', new_element % num_alphabet)
                print('Row ', i, ':\n', row)
        evolution_matrix.append(row)
        if debug == True:
            print('\nEvolution Matrix:\n', np.matrix(evolution_matrix))

    return evolution_matrix


def det_update_rule(num_elements, debug=False):
    """
    This function takes a string of numbers as input, which represent cells in a row. If a number exist outside the boundaries of the row, the string must be reentered.
        num_elements        -- Number of cells in row.
        valid               -- Flag that checks whether a string of numbers is valid.
    """

    valid = False

    while not valid:
        update_rule = []

        # Each number represents a cell in relation to the current one. +1 is one to the right. -1 is one to the left.
        string = input("Enter cells to be added for the rule (i.e. -2 0 3): ")
        update_rule = [int(d) for d in re.findall(r'-?\d+', string)]

        for i in update_rule:
            if abs(i) > num_elements:
                print(i,' is not a valid element for this automaton.\n')
                valid = False
                break
            else:
                valid = True     
        
        if debug == True:
            print('The update rule is ', update_rule)
    
    return update_rule


if __name__ == '__main__':
    """ Main function will serve as the driver for the program.
    It takes user input and runs cellular_automata() with given parameters:
            num_elements	-- Number of elements in automaton
            num_alphabet	-- Number of elements in the alphabet
                                -- In general, our alphabet will be numerical, starting at 0, and incrementing up to user input.
            start_state     -- Starting state of automaton.
			evolution_matrix     --n by n matrix which, when multiplied to a state in the automaton, will determine the next state.
			update_rule      -- Rule that governs how the automaton will evolve
                                    -- In general, the rule will be of the form: ( a + b + c + ... + n ) mod m, where a...n are elements of the previous state, and m is the number of elements in the alphabet.
            
            num_steps		-- Number of steps through which the automaton will evolve
    """

    # Debug mode -- displays math and logic of program where applicable
    debug = input('Debug? (Y/N) ')
    if debug == 'y' or debug == 'Y':
        debug = True
    else:
        debug = False

    num_elements = int(input('Enter # of elements in automaton: '))

    num_alphabet = int(input('Enter # of elements in alphabet: '))

    alphabet = []  # Used to make a nice visualization in the terminal

    for i in range(0, num_alphabet):
        alphabet.append(i)
    print('Your alphabet is', *alphabet)
    
    start_state = check_state(num_elements)

    update_rule = det_update_rule(num_elements, debug)

    evolution_matrix = evolve_matrix(num_elements, update_rule, debug)

    num_steps = int(input('Enter # of steps the automaton will take: '))

    print('\n\nBeginning process...\n')
    
    a = np.array(evolution_matrix)
    b = np.array([0,1,0,0])
    print(a, ' * ', b, ' = ', np.matmul(a,b))
    
    ### Currently using old non-evolution_matrix function
    #cellular_automata(num_elements, num_alphabet, start_state, num_steps, evolution_matrix)
