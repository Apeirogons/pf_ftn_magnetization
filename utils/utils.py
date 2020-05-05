
# Function to parse the other known sites found in other_known_sites.txt
def parse_other_known_sites(other_known_sites):
    splitted = [line.split(" ") for line in other_known_sites.split("\n")]
    known_sites = [(site[1],int(site[0])) for site in splitted]

    return known_sites

# Function to parse the mutations found in mutations.txt
def parse_mutations(mutations_text):
    un_spytagged = mutations_text.replace("SpyTag(N-term)","X00X")[:-1]
    m_removed = [line.split(" ")[1] for line in un_spytagged.split("\n")]
    plus_splitted = [line.split("+") for line in m_removed]
    parsed_mutations = []
    for mutations in plus_splitted:
        m = []
        for mutation in mutations:
            m.append((mutation[0], int(mutation[1:-1]), mutation[-1]))
        parsed_mutations.append(m)
    return parsed_mutations

# Function meant to find the aligned index of the original mutation site (as adding gaps in the alignment will shift sites to the right). 
# Currently unused, since there are no such mutations which require this
def get_true_site(site_n, seq_original):
    countdown = site_n-1
    true_site = -1
    for i, residue in enumerate(seq_original):
        if countdown == 0:
            true_site = i
            break
        if residue != "-":
            countdown -= 1
    return true_site   
