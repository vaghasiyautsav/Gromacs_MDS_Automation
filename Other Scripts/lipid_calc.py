def calculate_lipid_composition(cholesterol_percentage, lipid_percentage, total_lipids):
    cholesterol_count = int(total_lipids * (cholesterol_percentage / 100))
    lipid_count = total_lipids - cholesterol_count
    return cholesterol_count, lipid_count

def main():
    # List of possible lipids with details
    lipids = [
        {"name": "POPC", "area_per_lipid": 0.64, "details": "1-palmitoyl-2-oleoyl-glycero-3-phosphocholine. Model lipid for biophysical experiments."},
        {"name": "POPE", "area_per_lipid": 0.62, "details": "1-palmitoyl-2-oleoyl-glycero-3-phosphoethanolamine. Common in bacterial membranes."},
        {"name": "POPS", "area_per_lipid": 0.65, "details": "1-palmitoyl-2-oleoyl-glycero-3-phospho-L-serine. Negatively charged phospholipid."},
        {"name": "DOPC", "area_per_lipid": 0.72, "details": "1,2-dioleoyl-sn-glycero-3-phosphocholine. Used in artificial membranes."},
        {"name": "DOPE", "area_per_lipid": 0.63, "details": "1,2-dioleoyl-sn-glycero-3-phosphoethanolamine. Used in liposome formulation."}
    ]
    
    # Display the list of lipids for the user to select
    print("Select the type of lipid you want to use in the bilayer:")
    for i, lipid in enumerate(lipids, 1):
        print(f"{i}. {lipid['name']} - {lipid['details']}")
    
    # Ask user to select a lipid
    lipid_choice = int(input("Enter the number corresponding to your choice: "))
    if 1 <= lipid_choice <= len(lipids):
        selected_lipid = lipids[lipid_choice - 1]
    else:
        print("Invalid choice. Exiting.")
        return

    # Ask user for the composition of cholesterol and lipids
    cholesterol_percentage = float(input("Enter the percentage of cholesterol (e.g., 25 for 25%): "))
    lipid_percentage = 100 - cholesterol_percentage  # Assuming the rest is lipid

    # Ask user for the dimensions of the bilayer
    length = float(input("Enter the length of one side of the square bilayer (in Å): "))

    # Calculate the total number of lipids based on the area per lipid
    area_per_lipid = selected_lipid["area_per_lipid"]  # nm^2
    total_area = (length * length) / 100  # Convert Å^2 to nm^2
    total_lipids = int(total_area / area_per_lipid)

    # Calculate the number of cholesterol and lipid molecules
    cholesterol_count, lipid_count = calculate_lipid_composition(cholesterol_percentage, lipid_percentage, total_lipids)

    # Print the results
    print(f"\nFor a bilayer with dimensions {length} Å x {length} Å:")
    print(f"Total number of lipids: {total_lipids}")
    print(f"Number of cholesterol molecules: {cholesterol_count}")
    print(f"Number of {selected_lipid['name']} molecules: {lipid_count}")

if __name__ == "__main__":
    main()

