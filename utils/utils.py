

def parse_other_known_sites(other_known_sites):
    """
    Function to parse the other known sites found in other_known_sites.txt
    Args:
        other_known_sites: str (read from other_known_sites.txt)
    Returns:
        list of tuple(char, int): known_sites in the form of (residue, site number)
    """
    splitted = [line.split(" ") for line in other_known_sites.split("\n")]
    known_sites = [(site[1],int(site[0])) for site in splitted]

    return known_sites

def parse_mutations(mutations_text):
    """
    Function to parse the mutations found in mutations.txt
    Args:
        mutations_text: str (read from mutations.txt)
    Returns:
        list of tuple (char, int, char): mutations in the form of (original residue, site number, new residue)
    """
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

