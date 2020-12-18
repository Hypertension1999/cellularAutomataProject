import Nayuki
import numpy as np
import math
import sys
from hashlib import md5

"""*** Notes from meeting on 11/24: Don't need to double count cycles already counted. """
"""5mod5: Need to extend for longer cycle detection. Bardzell says there is another / more cycles that we have not shown."""
"""7mod3: Should output one cycle and one other cycle 20 < C < 90, but not all the other cycles. """
"""9mod3: Extreme irrevesible system that should only produce 1 1-cycle (after certain step, should evolve to 0-state forever)."""


def detect_cycle_transition(transition, alphabet, size):
    """
    Find T^k = T^n-k (IF there is one)
    finds the max number of cycle of transition matrix to compute
    """

    u_bound = 100000
    power = 1
    M_list = []
    ZERO = np.zeros((size, size), dtype=int)
    I = np.identity(size, dtype=int)

    result_matrix = transition
    for i in range(u_bound):
        M_list.append(result_matrix)
        result_matrix = (np.matmul(transition, result_matrix)) % alphabet
        power += 1

    for i in range(len(M_list)):
        for j in range(i + 1, len(M_list)):
            if i != j:
                # If true all states evloves to Zero State
                if(M_list[j] == I).all():
                    return(i, j, "I")
                if(M_list[j] == ZERO).all():
                    return(i, j, "0")
                if (M_list[i] == M_list[j]).all():
                    return(i, j, "Cycle")
            elif i == len(M_list):
                return(-1, "-1")
        else:
            continue

    return(-1, "-1")

def detect_cycle_transition_recurse_1(transition, alphabet, size):
    """
    Find T^k = T^n-k (IF there is one)
    finds the max number of cycle of transition matrix to compute
    """

    u_bound = 100
    M_List = []
    ZERO = np.zeros((size, size), dtype=int)
    I = np.identity(size, dtype=int)

    I_hash = md5(I.data)
    ZERO_hash = md5(ZERO.data)


    result_matrix = transition
    for i in range(u_bound):
        foo = md5(result_matrix.data)
        M_List.append(foo)
        result_matrix = (np.matmul(transition, result_matrix)) % alphabet

    #result_matrix = transition
    #M_List.append(result_matrix)
    return(detect_cycle_transition_recurse_2(M_List, transition, result_matrix, size, alphabet, ZERO_hash, I_hash, u_bound))

def detect_cycle_transition_recurse_2(M_List, transition, result_matrix, size, alphabet, ZERO_hash, I_hash, u_bound):

    for i in range(len(M_List)):
        for j in range(i + 1, len(M_List)):
            if i != j:
                # If true all states evloves to Zero State
                if(M_List[j].hexdigest() == I_hash.hexdigest()):
                    return(i, j, "I")
                if(M_List[j].hexdigest() == ZERO_hash.hexdigest()):
                    return(i, j, "0")
                if (M_List[i].hexdigest() == M_List[j].hexdigest()):
                    return(i, j, "Cycle")

    result_matrix = transition
    for i in range(u_bound):
        foo = md5(result_matrix.data)
        M_List.append(foo)
        result_matrix = (np.matmul(transition, result_matrix)) % alphabet

    #result_matrix = (np.matmul(transition, result_matrix)) % alphabet
    #M_List.append(result_matrix)
    return(detect_cycle_transition_recurse_2(M_List, transition, result_matrix, size, alphabet, ZERO_hash, I_hash, u_bound))


def is_reversible(B, rows, cols, size):
    """
    Find if N(T) == O-matrix
    """
    B.reduced_row_echelon_form()
    result = np.zeros([rows, cols], dtype=int)
    I = np.identity(size, dtype=int)

    for i in range(rows):
        for j in range(cols):
            result[i][j] = B.get(i, j)

    """
    if reduced_row_echelon_form of transition matrix = Identity matrix
    the matrix is reversible
    else it is irreversible
    """
    if (result == I).all():
    #     print("reversible")
        return(True)

    # print("irreversible")
    return(False)


def detect_unique_cycle(Nullspace_list, alphabet, size):
    """
        Find N(T^k -I) = cycles in automata
    """

    # taking an input list
    unique_nullspace = { "nullspace": [], "power": [] }

    # traversing the dict
    idx = 0
    for item in Nullspace_list["nullspace"]:
        if item not in unique_nullspace["nullspace"]:
            unique_nullspace["nullspace"].append(item)
            unique_nullspace["power"].append(Nullspace_list["power"][idx])
        idx += 1

    return(unique_nullspace)



def generate_automata_stats(data, size, alphabet):
    """
    RREF & NULL
    """

    data = data.tolist()
    rows = cols = size
    F = Nayuki.PrimeField(alphabet)
    B = Nayuki.Matrix(rows, cols, F)

    # Nayuki Matrix "B"
    for i in range(size):
        for j in range(size):
            B.set(i, j, data[i][j])

    # Transition Matrix "transition"
    transition = np.zeros([rows, cols], dtype=int)
    result_matrix = np.zeros([rows, cols], dtype=int)

    for i in range(rows):
        for j in range(cols):
            transition[i][j] = B.get(i, j)

    I = np.identity(size, dtype=int)        # Identity Matrix
    power = 1                               # power
    reversible = is_reversible(B, rows, cols, size)   # is automate reversible
    s, n, cycle_type = detect_cycle_transition_recurse_1(transition, alphabet, size)  # max steps

    #print("\nrref for matrix:")
    #B.reduced_row_echelon_form()
    #print(B)

    #print("\nnullspace for matrix:")
    #Basis = B.get_nullspace()
    #print(Basis)
    print("AAAAAAAAAAA")
    result_matrix_pow = transition
    Nullspace_list = { "nullspace": [], "power": [] }
    power = 1
    for i in range(n):
        #print("\n(T)^{} - I: ".format(power))
        #result_matrix = ( np.linalg.matrix_power( transition, power ) - I ) % alphabet

        if(i > 0):
            result_matrix_pow = (np.matmul(transition, result_matrix_pow)) % alphabet
        result_matrix = (result_matrix_pow - I) % alphabet

        # Set Nayuki Matrix to reult of (T)^n - I
        for i in range(size):  # For each column
            for j in range(size):
                B.set( i,j,int( result_matrix[i,j] ))

        #print("\nrref for (T)^{} - I: ".format(power))
        B.reduced_row_echelon_form()
        #print("\nnullspace for (T)^{} - I: ".format(power))
        Basis = B.get_nullspace()

        Nullspace_list["nullspace"].append(Basis)
        Nullspace_list["power"].append(power)

        power += 1

    # Get Automata Stats
    unique_nullspace = detect_unique_cycle( Nullspace_list, alphabet, size )  # max steps
    Automata_stats = []
    print("AAAAAAAAAAA")
    # Add final Cycle since the entire system is reversible
    if(reversible):
        result_matrix = I
        unique_nullspace["nullspace"].append(result_matrix)
        unique_nullspace["power"].append(power)

    # Prints all of the unique nullspaces to a file.
    for i in range( len(unique_nullspace["nullspace"]) ):
        power  = unique_nullspace["power"][i]
        cycles_size = len( unique_nullspace["nullspace"][i] )
        states = pow( alphabet, cycles_size )
        Automata_stats.append( { "nullspace": [], "power": 0, "cycles_size": 0, "cycles_count": 0, "states": 0} )

        subtract_repeat_states = 0
        if(i != 0):
            # Checks for divisibility
            for j in range(i, 0, -1):
                if j > 0 and power % unique_nullspace["power"][j-1] == 0:
                    item = Automata_stats[j-1]
                    subtract_repeat_states += item["states"]
            states -= subtract_repeat_states


        Automata_stats[i]["nullspace"] = unique_nullspace["nullspace"][i]
        Automata_stats[i]["power"] = power
        Automata_stats[i]["cycles_size"] = cycles_size
        Automata_stats[i]["cycles_count"] = states / power
        Automata_stats[i]["states"] = states

        if reversible:
            reversibility = 'Reversible system'
        else:
            reversibility = 'Irreversible system'

    return (Automata_stats, reversibility, s, n, cycle_type)
