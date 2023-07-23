def phytozomePacID(gff_file, ids_file):
    """
    _summary_
    this function analyzes the phytozome file for the 
    search of the corresponding pacid and provided a 
    agi_id file, phytozome file, it will search for the 
    corresponding pacid. 
    Arguments:
        gff_file -- _description_
        phytozome gff file and you can download from the phytozome 
        "Athaliana_447_Araport11.gene_exons.gff3"
        ids_file -- _description_
        takes the agiIDs for the search
    Returns:
        _description_
        returns a nested list with agi listed at 
        list[0] and the pacid listed at list[1].
    """
    phytozomedataframe = pd.read_csv(gff_file, sep = "\t")
    mRNA = phytozomedataframe.iloc[::,[2,8]]. \
                  where(phytozomedataframe.iloc[::,[2,8]]["gene"] == "mRNA").dropna()
    name = [i.split("=")[1] for i in ([j for i in ([i.split(";") \
                        for i in (mRNA.iloc[::,1].to_list())]) \
                                for j in i if j.startswith("Name=")])]
    pacid = [j for i in ([i.split(";") \
                          for i in (mRNA.iloc[::,1].to_list())]) \
                                    for j in i if j.startswith("pacid=")]
    agiPacID = [(i,j) for i,j in zip(name,pacid)]
    agi_ids = []
    final_ids = []
    with open(os.path.abspath(os.path.join(os.getcwd(), ids_file)), "r") as ids:
        for line in ids.readlines():
            agi_ids.append(line.strip())
        final_ids = [i.upper() for i in agi_ids if i!= ""]
        return [i for i in agiPacID for j in final_ids if j==i[0]]
