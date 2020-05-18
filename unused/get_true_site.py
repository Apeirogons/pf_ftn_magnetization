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
