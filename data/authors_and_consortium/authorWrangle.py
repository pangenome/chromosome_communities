import csv

with open("Authorship.csv", newline='') as f:
    reader = csv.DictReader(f)
    affiliation_dict, affiliation_dict_r = {}, {}
    affiliation_index = 1
    contribution_dict, contribution_dict_r = {}, {}
    contribution_index = 1
    authors = []
    authors_by_last_name = []
    funding = []
    
    authors_str = ''
    authors_by_last_name_str = ''
    for r in reader:
        last_name = r['NAME'].split()[-1].strip()
        first_name = " ".join(r['NAME'].split()[:-1]).strip()
        affiliation = r['AFFILIATION'].strip()
        contribution = r['CONTRIBUTION'].strip()
        rank = float(r["RANK"])
        rank = rank if rank != 10 else 3 # Make the 10s appear in the middle
        funding_string = r["FUNDING"]
        if funding_string != "":
            funding.append(funding_string.strip())
        if rank < 1:
            authors.append((rank, last_name, first_name, affiliation, contribution))
        if r["MEMBER"] == 'Yes':
            authors_by_last_name.append((last_name, first_name, rank, affiliation, contribution))

    authors.sort()
    HPRC_BANNER = False # True
    for rank, last_name, first_name, affiliations_string, contributions_string in authors:
        if rank == 10 and HPRC_BANNER:
            print("\nHuman Pangenome Reference Consortium")
            HPRC_BANNER = False

        # Work out affliations
        affiliations = [i.strip() for i in affiliations_string.split(";")]
        affiliation_indices = []
        for affiliation in affiliations:
            if affiliation not in affiliation_dict:
                affiliation_dict[affiliation] = affiliation_index
                affiliation_dict_r[affiliation_index] = affiliation
                affiliation_index += 1
            affiliation_indices.append(str(affiliation_dict[affiliation]))

        # Work out contributions
        contributions = [i.strip() for i in contributions_string.split(",")]
        for contribution in contributions:
            if len(contribution) > 0:
                if contribution not in contribution_dict:
                    contribution_dict[contribution] = []
                contribution_dict[contribution].append(first_name + " " + last_name)

        print("{} {}{},".format(first_name, last_name, ",".join(affiliation_indices)))
        authors_str += "{} {}{},".format(first_name, last_name, ",".join(affiliation_indices)) + ' '

    authors_by_last_name.sort()
    HPRC_BANNER = False # True
    for last_name, first_name, rank, affiliations_string, contributions_string in authors_by_last_name:
        if rank == 10 and HPRC_BANNER:
            print("\nHuman Pangenome Reference Consortium")
            HPRC_BANNER = False

        # Work out affliations
        affiliations = [i.strip() for i in affiliations_string.split(";")]
        affiliation_indices = []
        for affiliation in affiliations:
            if affiliation not in affiliation_dict:
                affiliation_dict[affiliation] = affiliation_index
                affiliation_dict_r[affiliation_index] = affiliation
                affiliation_index += 1
            affiliation_indices.append(str(affiliation_dict[affiliation]))

        # Work out contributions
        contributions = [i.strip() for i in contributions_string.split(",")]
        for contribution in contributions:
            if len(contribution) > 0:
                if contribution not in contribution_dict:
                    contribution_dict[contribution] = []
                contribution_dict[contribution].append(first_name + " " + last_name)

        print("{} {}{},".format(first_name, last_name, ",".join(affiliation_indices)))
        authors_by_last_name_str += "{} {}{},".format(first_name, last_name, ",".join(affiliation_indices)) + ' '


    for i in range(1, affiliation_index):
        print(i, affiliation_dict_r[i])

    #print("\nContributions")
    #for contribution in contribution_dict:
    #    print(contribution, "\n\t", ", ".join(contribution_dict[contribution]))

    #print("Funding\n", "\n".join(funding))
    
    print(authors_str)
    print('------------------------------')
    print(authors_by_last_name_str)
