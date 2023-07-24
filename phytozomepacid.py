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

#file version and two additional functions to prepare the files for the phytozome analysis
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
    with open(os.path.abspath(os.path.join(os.getcwd(),gff_file)), "r") as phytozome:
        with open(os.path.abspath(os.path.join(os.getcwd(),gff_file + "name")), "w") as phytozomer:
            for line in phytozome.readlines():
                if line.startswith("!"): 
                    continue     
            phytozomer.write(line)
    phytozomedataframe = pd.read_csv(os.path.abspath(os.path.join(os.getcwd(),gff_file + "name")), sep = "\t")
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

def prepareFunctionalNamePhytozome(phytozome_file, output_file):
    """
    _summary_
    this function prepares the files for the pacID analysis from the phytozome
    gene names files. Example file: Athaliana_447_Araport11.geneName.txt 
    Arguments:
        phytozome_file -- _description_
        phytozome file which contains the gene names and the description.
    Return: return_description    
        output_file -- _description_
        a nested list with the pacid and the gene names along with the splice variants. 
    """
    functionalName = {}
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file)), "r") as phytozome:
        with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file + "name")), "w") as phytozomer:
            for line in phytozome.readlines():
                if line.startswith("!"): 
                    continue     
            phytozomer.write(line)
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file + "name"), "r")) as file:
        for line in file.readlines():
            functionalName[line.strip().split("\t")[0]] = ''.join(j for i in \
                                                        ([line.strip().split("\t")[2:]]) for j in i)
    with open(os.path.abspath(os.path.join(os.getcwd(),output_file, "w") as processed:
        print(functionalName, file=processed)

def preparegeneNamePhytozome(phytozome_file, output_file):
    """
    _summary_
    this function prepares the files for the pacID analysis from the phytozome
    functional names files. Example file: Athaliana_447_Araport11.defline.txt 
    Arguments:
        phytozome_file -- _description_
        phytozome file which contains the gene names and the description.
    Return: return_description    
        output_file -- _description_
        a nested list with the pacid and the gene names along with the splice variants. 
    """   
    geneName = {}
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file)), "r") as phytozome:
        with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file + "name")), "w") as phytozomer:
            for line in phytozome.readlines():
                if line.startswith("!"): 
                    continue     
            phytozomer.write(line)
    with open(os.path.abspath(os.path.join(os.getcwd(),phytozome_file + "name")), "r")) as file:
        for line in file.readlines():
            geneName[line.strip().split("\t")[0]] = [line.strip().split("\t")[1]]
    with open(os.path.abspath(os.path.join(os.getcwd(),output_file), "w")) as processed:
        print(geneName, file=processed)
