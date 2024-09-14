def get_atom_name(amino_acid):
    # Define common atom names for salt bridge interactions
    atom_names = {
        'D': ['OD1', 'OD2'],  # Aspartic Acid
        'E': ['OE1', 'OE2'],  # Glutamic Acid
        'R': ['NH1', 'NH2'],  # Arginine
        'K': ['NZ'],          # Lysine
        'N': ['ND2'],         # Asparagine
        'Q': ['NE2'],         # Glutamine
        'H': ['ND1', 'NE2'],  # Histidine
    }
    return atom_names.get(amino_acid, ['CA'])  # Default to C-alpha if not found

def create_pairs():
    pairs = {}
    
    while True:
        print("Enter residue details (or type 'done' to finish):")
        
        res1 = input("First residue (e.g., D69): ")
        if res1.lower() == 'done':
            break
        res2 = input("Second residue (e.g., N300): ")
        
        # Extract amino acid type and number
        aa1, num1 = res1[0], res1[1:]
        aa2, num2 = res2[0], res2[1:]
        
        # Get atom names based on amino acid type
        atoms1 = get_atom_name(aa1)
        atoms2 = get_atom_name(aa2)
        
        # Create pairs
        for atom1 in atoms1:
            for atom2 in atoms2:
                pair_name = f"{num1}_{atom1}_{num2}_{atom2}"
                pair_value = f"resid {num1} and name {atom1} plus resid {num2} and name {atom2}"
                pairs[pair_name] = pair_value
                print(f"Pair added: {pair_name} = \"{pair_value}\"")
    
    print("\nAll pairs created:")
    for key, value in pairs.items():
        print(f"pairs[\"{key}\"]=\"{value}\"")

create_pairs()
