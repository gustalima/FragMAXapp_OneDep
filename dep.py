import os
from glob import glob
import sys
import shutil
import subprocess


path = os.getcwd()
group_dep_id = "G_1002167"

protein = "PrtK"
protein_full_name = "Proteinase K"
library_full_name = "Frag Xtal Screen"
library = "JBS"

taxon_ID = "37998"

keywords = "fragment screening, hydrolase, inhibition"
protein_function = "Hydrolase"

citation_title = "To be published"
journal_abrev = "To be published"
journal_vol = "?"
first_page = "?"
last_page = "?"
year = "?"
PUBMED_id = "?"
DOI = "?"

collection_date = "2019-10-22"

organism_scientific = "'Parengyodontium album'"
nat_strain = "?"

gene_src_gene = "PROK"

host_org_scientific_name = "'Escherichia coli'"
host_org_ncbi_taxonomy_id = "562"

database_accession = "P06873"
database_accession_uniprot = "P06873"

beg_seq_num = "1"
end_seq_num = "279"

source_method = "man"

sequence = "AAQTNAPWGLARISSTSPGTSTYYYDESAGQGSCVYVIDTGIEASHPEFEGRAQMVKTYYYSSRDGNGHGTHCAGTVGSRTYGVAKKTQLFGVKVLDDNGSGQYSTIIAGMDFVASDKNNRNCPKGVVASLSLGGGYSSSVNSAAARLQSSGVMVAVAAGNNNADARNYSPASEPSVCTVGASDRYDRRSSFSNYGSVLDIFGPGTDILSTWIGGSTRSISGTSMATPHVAGLAAYLMTLGKTTAASACRYIADTANKGDLSNIPFGTVNLLAYNNYQA"

crystallisation_condition = "1.2M ammonium sulfate, 0.1M Tris-HCl, pH 8.0"
crystallisation_temp = "290"
crystallisation_ph = "8.0"
crystallisation_method = "VAPOR DIFFUSION, SITTING DROP"

EC = "3.4.21.64"

### DEFINITIONS FOR NEW ENTRIES

u = {
    1: {
        "authorID": "1",
        "salutation": "Dr.",
        "firstName": "Gustavo",
        "lastName": "Lima",
        "middleName": "'M A'",
        "role": "'principal investigator/group leader'",
        "orgType": "academic",
        "email": "gustavo.lima@maxiv.lu.se",
        "address": "'Fotongatan 2'",
        "city": "Lund",
        "state": "Skane",
        "postalCode": "22484",
        "Country": "Sweden",
        "faxNumber": "?",
        "phoneNumber": "'+46 725 24 6006'",
        "ORCID": "0000-0002-0011-5559",
    },
    2: {
        "authorID": "2",
        "salutation": "Dr.",
        "firstName": "Vladimir",
        "lastName": "Talibov",
        "middleName": "?",
        "role": "'responsible scientist'",
        "orgType": "academic",
        "email": "vladimir.talibov@maxiv.lu.se",
        "address": "'Fotongatan 2'",
        "city": "Lund",
        "state": "Skane",
        "postalCode": "22484",
        "Country": "Sweden",
        "faxNumber": "?",
        "phoneNumber": "'+00 000 00 0000'",
        "ORCID": "0000-0002-1135-2744",
    },
    3: {
        "authorID": "3",
        "salutation": "Ms.",
        "firstName": "Laila",
        "lastName": "Benz",
        "middleName": "S",
        "role": "'responsible scientist'",
        "orgType": "academic",
        "email": "laila.benz@fu-berlin.de",
        "address": "'Institut fur Biochemie und Chemie Thielallee 63'",
        "city": "Berlin",
        "state": "Berlin",
        "postalCode": "14195",
        "Country": "Germany",
        "faxNumber": "?",
        "phoneNumber": "'+00 000 00 0000'",
        "ORCID": "0000-0001-6369-8707",
    },
    4: {
        "authorID": "4",
        "salutation": "Mr.",
        "firstName": "Elmir",
        "lastName": "Jagudin",
        "middleName": "?",
        "role": "investigator",
        "orgType": "academic",
        "email": "elmir.jagudin@maxiv.lu.se",
        "address": "'Fotongatan 2'",
        "city": "Lund",
        "state": "Skane",
        "postalCode": "22484",
        "Country": "Sweden",
        "faxNumber": "?",
        "phoneNumber": "'+00 000 00 0000'",
        "ORCID": "0000-0001-7005-7641",
    },
    5: {
        "authorID": "5",
        "salutation": "Dr.",
        "firstName": "Uwe",
        "lastName": "Mueller",
        "middleName": "?",
        "role": "investigator",
        "orgType": "academic",
        "email": "uwe.muller@maxiv.lu.se",
        "address": "'Fotongatan 2'",
        "city": "Lund",
        "state": "Skane",
        "postalCode": "22484",
        "Country": "Sweden",
        "faxNumber": "?",
        "phoneNumber": "'+00 000 00 0000'",
        "ORCID": "0000-0002-7139-0718",
    },
}


_software = """# 
loop_
_software.pdbx_ordinal 
_software.name 
_software.version 
_software.date 
_software.type 
_software.contact_author 
_software.contact_author_email 
_software.classification 
_software.location 
_software.language 
1 XSCALE      ?        ?              package 'Wolfgang Kabsch'    ?                        
'data scaling'    http://www.mpimf-heidelberg.mpg.de/~kabsch/xds/html_doc/xscale_program.html ? 
2 REFMAC      5.8.0238 ?              program 'Garib N. Murshudov' garib@ysbl.york.ac.uk    
refinement        http://www.ccp4.ac.uk/dist/html/refmac5.html        Fortran_77 
3 PDB_EXTRACT 3.25     'Apr. 1, 2019' package PDB                  deposit@deposit.rcsb.org 
'data extraction' http://sw-tools.pdb.org/apps/PDB_EXTRACT/           C++        
4 XDS         .        ?               program ?                    ?                        
'data reduction'  ? ?          
5 PHENIX      1.16.3549 ?              program ? ?    
refinement        ?        ?    
6 REFMAC      5.8.0238 ?              program 'Garib N. Murshudov' garib@ysbl.york.ac.uk    
phasing        http://www.ccp4.ac.uk/dist/html/refmac5.html        Fortran_77   
"""


#######


def generate_data_template(pdbFile, fragment, sub, run):
    path = "/".join(pdbFile.split("/")[:-1])
    protein = pdbFile.split("/")[-1].split("-")[0]
    library = pdbFile.split("/")[-1].split("-")[1]
    pdb_model = f"{protein}-{library}-{fragment}{sub}_{run}.pdb"
    print(pdb_model)
    with open(f"{path}/{pdb_model}", "r") as readFile:
        pdb = readFile.readlines()
    for line in pdb:
        if "CRYST1" in line:
            space_group = " ".join(line.split()[7:])
            a, b, c, alpha, beta, gamma = line.split()[1:7]
        if "PHENIX" in line:
            sw = "PHENIX"
        if "REFMAC" in line:
            sw = "REFMAC"
    data_template = f"""

<contact_author_PI_id= 1 >           !(must be given 1)
<contact_author_PI_salutation= Dr. >    !(Dr./Prof./Mr./Mrs./Ms.)
<contact_author_PI_first_name= John >    !(e.g. John)
<contact_author_PI_last_name= Rodgers >     !(e.g. Rodgers)
<contact_author_PI_middle_name= X >         
<contact_author_PI_role= principal investigator/group leader> !Fixed. Do not change!
<contact_author_PI_organization_type= academic>  !(or commercial, government, other)
<contact_author_PI_email= name@host.domain.country >              !(e.g.   name@host.domain.country)      
<contact_author_PI_address= 610 Taylor road >            !(e.g.  610 Taylor road)
<contact_author_PI_city= Piscataway >               !(e.g.  Piscataway)
<contact_author_PI_State_or_Province= New Jersey >  !(e.g.  New Jersey)
<contact_author_PI_Zip_Code= 08864 >           !(e.g.  08864)
<contact_author_PI_Country= United States >        !(e.g. United States, United Kindom, . )
<contact_author_PI_fax_number=  >
<contact_author_PI_phone_number= +00 000 00 0000 >    !(e.g.  01(617) 555-1213 )


<Release_status_for_coordinates= HOLD FOR PUBLICATION >      !(e.g. HOLD FOR PUBLICATION)
<Release_status_for_structure_factor= HOLD FOR PUBLICATION > !(e.g. HOLD FOR PUBLICATION)
<Release_status_for_sequence= HOLD FOR RELEASE >             !(RELEASE NOW  or  HOLD FOR RELEASE)

<structure_title= Structure_title >     !(e.g. Crystal Structure Analysis of the B-DNA)
<structure_details= Structure_details >  

<structure_author_name= Rodgers, J. >

<primary_citation_author_name= Rodgers, J. >


<primary_citation_id= primary>     
<primary_citation_journal_abbrev= title of publication or manuscript >     (e.g. To be published)
<primary_citation_title= To be published >   
<primary_citation_year= ? >
<primary_citation_journal_volume= ? > 
<primary_citation_page_first= ? >
<primary_citation_page_last= ? >

<molecule_id= 1 >        (e.g. 1 )
<molecule_name= {protein} >       (e.g.  RNA Hammerhead Ribozyme )
<molecule_type= polymer >    (e.g. polymer , non-polymer, macrolide  )
<molecule_source_method= {source_method} >   (e.g. man , nat, syn)

<Molecular_entity_id= 1 >       (e.g. 1 )
<Enzyme_Comission_number= {EC} >   (if known: e.g. 2.7.7.7)

<natural_source_entity_id= 1 >          (e.g. 1, 2..)
<natural_source_scientific_name= {organism_scientific} >    (e.g. Homo sapiens)
<natural_source_organism_strain= {nat_strain} >    (e.g. DH5a , BMH 71-18)
<natural_source_details= cell >            (e.g. organ, tissue, cell ..)

<structure_keywords= {protein_function} >  !(e.g. beta barrel)

<biological_assembly_chain_number= 1 >  !(e.g.  1 for monomer, 2 for dimer ..)

<crystal_number= 1 >            (e.g. 1, )
<crystallization_method= {crystallisation_method} >      !(e.g. BATCH MODE, EVAPORATION, SLOW COOLING) 
<crystallization_pH= {crystallisation_ph} >          (e.g. 7.5 ...)
<crystallization_temperature= {crystallisation_temp} > !(e.g. 298) (in Kelvin) 
<crystallization_details= {crystallisation_condition} >     !(e.g. 5% DMSO, 100 mM HEPES;  PEG 4000, NaCl etc.)

<radiation_experiment= 1 >      !(e.g. 1, 2, ...)
<radiation_source= SYNCHROTRON >           !(e.g. SYNCHROTRON, ROTATING ANODE ..)
<radiation_source_type= MAX IV BEAMLINE BIOMAX >      !(e.g. NSLS BEAMLINE X8C ..)
<radiation_wavelength_id= 1 >
<radiation_wavelengths= 0.827 >       !(e.g. 1.502, or a list 0.987,0.988 ..)
<radiation_detector= PIXEL  >         !(e.g. CCD, PIXEL, AREA DETECTOR, IMAGE PLATE ..)
<radiation_detector_type= DECTRIS EIGER X 16M >     !(e.g. CCD,  ADSC QUANTUM 1,  ..)
<radiation_detector_details=  >    (e.g. mirrors...)
<data_collection_date= {collection_date} >             !(e.g. 2004-01-07)
<data_collection_temperature= 100 >      !(e.g. 100 ) (in Kelvin)
<data_collection_protocol= SINGLE WAVELENGTH >          !(e.g. SINGLE WAVELENGTH, MAD, ...)
<data_collection_monochromator= Si(111) >     (e.g. GRAPHITE, Ni FILTER ...)
<data_collection_monochromatic_or_laue=  M >  !(default M, give L if Laue diffr.)



<refinement_detail=   >  
<refinement_start_model=   >    (e.g. pdbid 100D)

<database_entity_id= 1  >  (e.g. 1 )
<database_name= UNP >  (e.g. BMCD, BMRB, EMDB, PDB, NDB, TargetTrack )
<database_code= {database_accession_uniprot} >  (e.g. 1ABC, 100D, TNKS2_HUMAN )
<database_accession= {database_accession_uniprot}  >  (e.g. 100D, Q9H2K2  )

<space_group = {space_group}> (use International Table conventions)
<unit_cell_a     =  {a}  >
<unit_cell_b     =  {b}  >
<unit_cell_c     =  {c}  >
<unit_cell_alpha =  {alpha} >
<unit_cell_beta  =  {beta}  >
<unit_cell_gamma =  {gamma} >

<molecule_entity_id= 1 >
<molecule_entity_type= polypeptide(L) >
<molecule_one_letter_sequence= 
{sequence}
>
< molecule_chain_id=A >
< target_DB_id=  > (if known) 
< sequence_database_id= {database_accession_uniprot} > (if known) 
< sequence_database_name= UNP > (if known) 

"""

    with open(f"{path}/data_template_{fragment}{sub}_{run}.text", "w") as writeFile:
        writeFile.write(data_template)


def generate_log_script(pdbFile, fragment, sub, run):
    path = "/".join(pdbFile.split("/")[:-1])
    protein = pdbFile.split("/")[-1].split("-")[0]
    library = pdbFile.split("/")[-1].split("-")[1]
    pdb_model = f"{protein}-{library}-{fragment}{sub}_{run}.pdb"
    refine_log = f"{protein}-{library}-{fragment}{sub}_{run}.log"
    unmerged_data = f"{protein}-{library}-{fragment}{sub}_{run}_data_F.mtz"
    refined_data = f"{protein}-{library}-{fragment}{sub}_{run}.mtz"
    scale_sw_log = "CORRECT.LP"
    integrate_sw_log = "INTEGRATE.LP"

    with open(f"{path}/{pdb_model}", "r") as readFile:
        pdb = readFile.readlines()

    for line in pdb:
        if "PHENIX" in line:
            sw = "PHENIX"
        if "REFMAC" in line:
            sw = "REFMAC"

    with open(f"{path}/{pdb_model}", "w") as writeFile:
        new_pdb = [x for x in pdb if not x.startswith("TITLE") and not x.startswith("COMPND")]
        writeFile.write("".join(new_pdb))
    log_inp = f"""


<data_indexing_software = "XDS" >
<data_indexing_LOG_file_name = "{path}/{protein}{library}-{fragment}{sub}_{run}_{integrate_sw_log}" >
<data_indexing_CIF_file_name = " " >  (if mmCIF format)

<refine_software= {sw}  > (e.g. remfac, phenix)
<refine_xyz_file= {path}/{pdb_model}  > (coordinate file in PDB or mmcif format)
<refine_log_file= {path}/{refine_log}  > (optional: log file in PDB/mmcif format)

<data_scaling_software = "XSCALE" >
<data_scaling_LOG_file_name = "{path}/{protein}{library}-{fragment}{sub}_{run}_{scale_sw_log}" >
<data_scaling_CIF_file_name = " " >  (if  mmCIF format)

<data_template_file= {path}/data_template_{fragment}{sub}_{run}.text  > 

<mr_software= PHASER >  (e.g. PHASER, MOLREP ..)
<mr_log_file_LOG_1=   >  (log file containing statistics)
<mr_log_file_LOG_2=   >  (log file containing statistics)

<phasing_method = "FOURIER SYNTHESIS" >        
<phasing_software = " " >

<reflection_data_file_name= {path}/{refined_data}  >  (give SF file name)
<reflection_data_detail=  >  (optional, give a note to the data set)
<reflection_data_free_set= 1  >  (Free set: input a number if not 0)

<reflection_data_file_name= {path}/{unmerged_data}  >   (give SF file name)
<reflection_data_detail=   >  (optional, give a note to the data set)

<statistics_output= {protein}-{library}-{fragment}{sub}_{run}.mmcif  >    (for coordinates and statistics)
<sf_output= {protein}-{library}-{fragment}{sub}_{run}_sf.mmcif  >            (for structure factors)
    """
    with open(f"{path}/log_script_{fragment}{sub}_{run}.inp", "w") as writeFile:
        writeFile.write(log_inp)


def create_pdbx_contact_author(authorsDict):
    _pdbx_contact_author_HEADER = f"""# 
loop_
_pdbx_contact_author.id 
_pdbx_contact_author.name_salutation 
_pdbx_contact_author.name_first 
_pdbx_contact_author.name_last 
_pdbx_contact_author.name_mi 
_pdbx_contact_author.role 
_pdbx_contact_author.organization_type 
_pdbx_contact_author.email 
_pdbx_contact_author.address_1 
_pdbx_contact_author.city 
_pdbx_contact_author.state_province 
_pdbx_contact_author.postal_code 
_pdbx_contact_author.country 
_pdbx_contact_author.fax 
_pdbx_contact_author.phone 
_pdbx_contact_author.identifier_ORCID
"""

    _pdbx_contact_author_BODY = ""
    for k, authDict in authorsDict.items():
        entry_auth_body = f"""{authDict["authorID"]} {authDict["salutation"]} {authDict["firstName"]}   {authDict["lastName"]}       {authDict["middleName"]}     {authDict["role"]} {authDict["orgType"]} 
{authDict["email"]}         {authDict["address"]} {authDict["city"]}  {authDict["state"]} {authDict["postalCode"]} {authDict["Country"]} ? 
{authDict["phoneNumber"]} {authDict["ORCID"]} """
        _pdbx_contact_author_BODY += entry_auth_body
        _pdbx_contact_author_BODY += "\n"
    return _pdbx_contact_author_HEADER + _pdbx_contact_author_BODY


def parse_citation_author(authorsDict):
    citation_author_HEAD = """# 
loop_
_citation_author.citation_id 
_citation_author.name 
_citation_author.ordinal 
_citation_author.identifier_ORCID \n"""
    citation_author_BODY = ""
    for k, v in u.items():
        fN = v["firstName"][0]
        lN = v["lastName"]
        mN = "".join([x + "." for x in v["middleName"].split() if "?" not in x])
        mN = mN.replace("'", "").replace('"', "").replace("?", "")
        author = f"'{lN}, {fN}.{mN}'"
        largeName = max([(len(v["lastName"]) + 5 + len(v["middleName"].split())) for v in u.values()])
        blanks = largeName - len(author)
        citation_author_BODY += "primary " + author + (blanks + 1) * " " + v["authorID"] + " " + v["ORCID"] + "\n"

    return citation_author_HEAD + citation_author_BODY


def parse_audit_author(authorsDict):
    _audit_author_HEAD = """# 
loop_
_audit_author.pdbx_ordinal 
_audit_author.name  \n"""
    _audit_author_BODY = ""
    for k, v in u.items():
        fN = v["firstName"][0]
        lN = v["lastName"]
        mN = "".join([x + "." for x in v["middleName"].split() if "?" not in x])
        mN = mN.replace("'", "").replace('"', "").replace("?", "")
        author = f"'{lN}, {fN}.{mN}'"
        _audit_author_BODY += v["authorID"] + " " + author + "\n"

    return _audit_author_HEAD + _audit_author_BODY


def fix_occupancies_conformer(mmcif):
    with open(mmcif, "r") as r:
        _mmcif = r.readlines()
    coordinates = [x for x in _mmcif if x.startswith("ATOM") or x.startswith("HETATM")]
    atom_c = [x for x in coordinates if x.startswith("ATOM")]
    hetatm_c = [x for x in coordinates if x.startswith("HETATM")]

    res = atom_c[:20][0].split()[5]
    rDescription = dict()
    occupDict = dict()
    resD = dict()
    for atom_line in atom_c:
        values = atom_line.split()
        atomn = values[2]
        atomi = values[1]
        conf = values[4]
        resn = values[5]
        resi = values[10]
        occ = float(values[-6])
        resD[atomn + "_" + atomi + "_" + conf] = occ
        occupDict[resn + "_" + resi] = resD
        if resn + "_" + resi not in rDescription and len(rDescription) > 0:
            rDescription[resn + "_" + resi] = resD
            resD = dict()
        if resn + resi not in rDescription or rDescription[resn + "_" + resi] != "":
            rDescription[resn + "_" + resi] = ""

    need_fix = False
    conf_list = list()
    for k, v in occupDict.items():
        for k2, v2 in v.items():
            if "_A" in k2:
                # print(k,k2)
                # print([x for x in v if k2.replace("_A","") in x])
                atms = [x for x in v if k2.split("_")[0] == x.split("_")[0] and "." not in x]
                occ_conf = sum([v[x] for x in atms])
                if occ_conf > 1:
                    need_fix = True
                    # print(k,k2,atms,round(abs(occ_conf-1),2))
                    conf_list.append(atms + [round(abs(occ_conf - 1), 2)])

    if need_fix:

        foic = conf_list
        if foic is not None:
            for at in foic:
                for n, line in enumerate(_mmcif):
                    if line.startswith("ATOM"):
                        if all(ext in line for ext in at[-2].split("_")):
                            old_oc = f" {line.split()[-6]} "
                            new_oc = " {0:.2f} ".format(float(float(line.split()[-6]) - at[-1]))
                            print(f"old occ {old_oc}")
                            print(f"new occ {new_oc}")
                            _mmcif[n] = _mmcif[n].replace(old_oc, new_oc)
    return _mmcif


def fix_occupancies_ensemmble(mmcif):
    with open(mmcif, "r") as readFile:
        xml = readFile.readlines()
    occupDict = dict()
    residDict = dict()
    for line in xml:
        resid = ""
        if line.startswith("ATOM") or line.startswith("HETATM"):
            lineSplit = line.split()
            resid = lineSplit[1] + "_" + lineSplit[2] + "_" + lineSplit[5] + "_" + lineSplit[10] + "_" + lineSplit[8]
            resid2 = lineSplit[2] + "_" + lineSplit[5] + "_" + lineSplit[10] + "_" + lineSplit[8]
            occup = lineSplit[16]
            if float(occup) != 1.0:
                if resid in occupDict:
                    addedOccup = float(occupDict[resid]) + float(occup)
                    occupDict[resid] = addedOccup
                else:
                    occupDict[resid] = float(occup)
                if resid2 in residDict:
                    addedOccup = float(residDict[resid2]) + float(occup)
                    residDict[resid2] = addedOccup
                else:
                    residDict[resid2] = float(occup)

    for k, v in residDict.items():
        if v > 1.0:
            atomsInRes = [x for x in occupDict if k in x]
            tochange = atomsInRes[-1]
            changeOccup = round(v - 1.00, 2)
            if type(occupDict[tochange]) == float and type(changeOccup) == float:
                print("Fixing occupancy for " + tochange)
                new_occ = occupDict[tochange] - changeOccup
                if new_occ < 0:
                    new_occ = 0
                print(f"old occ {occupDict[tochange]}")
                print("new occ " + "{0:.2f}".format(new_occ))
                occupDict[tochange] = "{0:.2f}".format(new_occ)
    output_mmcif = ""
    for line in xml:
        resid = ""
        if line.startswith("ATOM") or line.startswith("HETATM"):

            lineSplit = line.split()
            resid = lineSplit[1] + "_" + lineSplit[2] + "_" + lineSplit[5] + "_" + lineSplit[10] + "_" + lineSplit[8]
            occup = lineSplit[16]
            if resid in occupDict:
                output_mmcif += line.replace(occup, f"{occupDict[resid]}")
                if "\n" not in line.replace(occup, f"{occupDict[resid]}"):
                    output_mmcif += "\n"
            else:
                output_mmcif += line
                if "\n" not in line:
                    output_mmcif += "\n"
        else:
            output_mmcif += line
            if "\n" not in line:
                output_mmcif += "\n"
    return output_mmcif


def generate_refine_ls_shell_table(final_pdb, logFile):
    output_table = list()
    output_table = [
        "#",
        "loop_",
        "_refine_ls_shell.d_res_high",
        "_refine_ls_shell.d_res_low",
        "_refine_ls_shell.pdbx_total_number_of_bins_used",
        "_refine_ls_shell.percent_reflns_obs",
        "_refine_ls_shell.number_reflns_R_work",
        "_refine_ls_shell.R_factor_all",
        "_refine_ls_shell.R_factor_R_work",
        "_refine_ls_shell.R_factor_R_free",
        "_refine_ls_shell.percent_reflns_R_free",
        "_refine_ls_shell.number_reflns_R_free",
        "_refine_ls_shell.R_factor_R_free_error",
        "_refine_ls_shell.number_reflns_all",
        "_refine_ls_shell.number_reflns_obs",
        "_refine_ls_shell.pdbx_refine_id",
    ]
    search_file = final_pdb
    if not select_search_file(search_file):
        search_file = logFile

    with open(search_file, "r") as readFile:
        file_content = readFile.readlines()
    if ".pdb" in search_file:
        file_format = "pdb"
    if ".log" in search_file:
        file_format = "reflog"

    if "refmac" in "".join(file_content).lower():
        table_ref = find_refine_table_refmac(file_content, file_format, logFile)
    else:
        table_ref = find_refine_table_phenix(file_content, file_format, logFile)

    for l in table_ref:
        completeness = "{0:.4f}".format(float(l[4]) * 100)
        binSize = str(len(table_ref))
        l = [str(x) for x in l]
        table_line = (
            l[3]
            + " "
            + l[1]
            + " "
            + binSize
            + " "
            + completeness
            + (9 - len(completeness)) * " "
            + l[5]
            + " "
            + "?"
            + " "
            + l[7]
            + " "
            + l[8]
            + " "
            + "?"
            + " "
            + l[6]
            + " "
            + "0.0000"
            + " "
            + str(int(l[5]) + int(l[6]))
            + "  ? 'X-RAY DIFFRACTION'"
        )

        output_table.append(table_line)
    return [x + "\n" for x in output_table]


def parse_statistics_from_logs(logFile, correctLP):
    statsDict = dict()

    with open(logFile, "r") as readFile:
        ref_log = readFile.readlines()
    if "refmac" in "".join(ref_log).lower():
        with open(logFile.replace(".log", ".pdb"), "r") as readFile:
            ref_pdb = readFile.readlines()
        for line in ref_pdb:
            if "R VALUE     (WORKING + TEST SET)" in line:
                ls_R_factor_R_free = line.split()[-1]
            if "R VALUE            (WORKING SET)" in line:
                ls_R_factor_R_work = line.split()[-1]
            if "NUMBER OF REFLECTIONS" in line:
                number_obs = line.split()[-1]
            if "FREE R VALUE TEST SET COUNT" in line:
                number_reflns_R_free = line.split()[-1]
                if number_reflns_R_free == "NULL":
                    number_reflns_R_free = 0
            if "RESOLUTION RANGE HIGH (ANGSTROMS)" in line:
                d_resolution_high = line.split()[-1]
    else:
        for line in ref_log:
            if "Final R-work" in line:
                ls_R_factor_R_free = line.split()[-1]
                ls_R_factor_R_work = line.split()[3]
            if "remove outliers" in line:
                number_obs = line.split()[-1]
            if "resolution" in line and "A, n_refl" in line and "% free" in line:
                rfree_perc = line.split()[-3]
                if "number_obs" in locals():
                    number_reflns_R_free = int(float(number_obs) * float(rfree_perc) / 100)
                d_resolution_high = line.split()[1]

    with open(correctLP, "r") as readFile:
        cor_log = readFile.readlines()
    for line in cor_log:
        if "    total" in line:
            percent_possible_obs = line.split()[4].replace("%", "")
            percent_possible_obs = "{0:.3f}".format(float(percent_possible_obs))
            pdbx_netI_over_sigmaI = "{0:.3f}".format(float(line.split()[8]))
            l = line.split()
            pdbx_CC_half = line.split()[-4].replace("*", "")
            pdbx_CC_half = "{0:.3f}".format(float(pdbx_CC_half) / 100)
            pdbx_Rmerge_I_obs = p2f(l[5])
            pdbx_Rrim_I_all = p2f(l[9])
            pdbx_number_measured_all = l[1]
        if "WILSON LINE (using all data)" in line:
            B_iso_Wilson_estimate = line.split()[-3]
        if "CHI^2-VALUE OF FIT OF CORRECTION FACTORS" in line:
            pdbx_chi_squared = line.split()[-1]
        if " NUMBER OF ACCEPTED OBSERVATIONS" in line:
            unique = line.split()[-1]
        if " NUMBER OF UNIQUE ACCEPTED REFLECTIONS" in line:
            accept = line.split()[-1]
    pdbx_redundancy = "{0:.3f}".format(float(unique) / float(accept))

    statsDict = {
        "_refine.ls_R_factor_R_free": ls_R_factor_R_free,
        "_refine.ls_R_factor_R_work": ls_R_factor_R_work,
        "_refine.ls_number_reflns_R_free": number_reflns_R_free,
        "_reflns.number_obs": number_obs,
        "_reflns.d_resolution_high": d_resolution_high,
        "_reflns.percent_possible_obs": percent_possible_obs,
        "_reflns.pdbx_netI_over_sigmaI": pdbx_netI_over_sigmaI,
        "_reflns.pdbx_Rmerge_I_obs": pdbx_Rmerge_I_obs,
        "_reflns.B_iso_Wilson_estimate": B_iso_Wilson_estimate,
        "_reflns.pdbx_redundancy": pdbx_redundancy,
        "_reflns.pdbx_Rrim_I_all": pdbx_Rrim_I_all,
        "_reflns.pdbx_CC_half": pdbx_CC_half,
        "_reflns.pdbx_number_measured_all": pdbx_number_measured_all,
        "_reflns.pdbx_chi_squared": pdbx_chi_squared,
    }
    return statsDict


def create_group_description(protein, library, group_dep_id):

    group_description = f"""# 
_pdbx_deposit_group.group_id            {group_dep_id} 
_pdbx_deposit_group.group_description 
;PanDDA analysis of {library_full_name} vs. {protein}, including auto-refined models with ligands placed according to PanDDA-map and automatically refined models necessary to reproduce ground state model 
; 
_pdbx_deposit_group.group_title 
'PanDDA analysis of {library_full_name} vs. {protein}' 
_pdbx_deposit_group.group_type          undefined """

    return group_description

    """with open(mmcif,"r") as readFile:
        a=readFile.read()
    if "#" in a[-3:]:
        group_description=group_description[1:]
    if "_pdbx_deposit_group.group_title" not in a:
        with open(mmcif,"w") as writeFile:
            new_mmcif=a+group_description
            writeFile.write(new_mmcif)"""


def parse_pdbx_entity_instance_feature():
    _pdbx_entity_instance_feature = """#
_pdbx_entity_instance_feature.ordinal        1 
_pdbx_entity_instance_feature.comp_id        XXX 
_pdbx_entity_instance_feature.asym_id        ? 
_pdbx_entity_instance_feature.seq_num        ? 
_pdbx_entity_instance_feature.auth_comp_id   XXX 
_pdbx_entity_instance_feature.auth_asym_id   ? 
_pdbx_entity_instance_feature.auth_seq_num   ? 
_pdbx_entity_instance_feature.feature_type   'SUBJECT OF INVESTIGATION' 
_pdbx_entity_instance_feature.details        ? """
    return _pdbx_entity_instance_feature


def retrive_reflns_shell_section(correctLP):
    n_line = (
        "loop_ \n"
        + "_reflns_shell.pdbx_diffrn_id \n"
        + "_reflns_shell.pdbx_ordinal \n"
        + "_reflns_shell.d_res_high \n"
        + "_reflns_shell.d_res_low \n"
        + "_reflns_shell.number_measured_obs \n"
        + "_reflns_shell.number_measured_all \n"
        + "_reflns_shell.number_unique_obs \n"
        + "_reflns_shell.pdbx_rejects \n"
        + "_reflns_shell.Rmerge_I_obs \n"
        + "_reflns_shell.meanI_over_sigI_obs \n"
        + "_reflns_shell.pdbx_Rsym_value \n"
        + "_reflns_shell.pdbx_chi_squared \n"
        + "_reflns_shell.pdbx_redundancy \n"
        + "_reflns_shell.percent_possible_obs \n"
        + "_reflns_shell.pdbx_netI_over_sigmaI_obs \n"
        + "_reflns_shell.number_possible \n"
        + "_reflns_shell.number_unique_all \n"
        + "_reflns_shell.Rmerge_F_all \n"
        + "_reflns_shell.Rmerge_F_obs \n"
        + "_reflns_shell.Rmerge_I_all \n"
        + "_reflns_shell.meanI_over_sigI_all \n"
        + "_reflns_shell.percent_possible_all \n"
        + "_reflns_shell.pdbx_Rrim_I_all \n"
        + "_reflns_shell.pdbx_Rpim_I_all \n"
        + "_reflns_shell.pdbx_CC_half \n"
        + "_reflns_shell.pdbx_number_anomalous \n"
    )

    with open(correctLP, "r") as readFile:
        cor = readFile.readlines()

    for n, line in enumerate(cor):
        if "SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE" in line:
            line1 = n
        if "total   " in line:
            linelast = n
        if "INCLUDE_RESOLUTION_RANGE" in line:
            d_res_low = line.split()[1]
            d_res_low = "{0:.2f}".format(float(d_res_low))

    l = cor[line1 + 4 : linelast]
    output = list()
    for n, j in enumerate(l, 1):
        if n != 1:
            d_res_low = l[n - 2].split()[0]

        k = j.split()

        pdbx_diffrn_id = 1
        pdbx_ordinal = n
        d_res_high = k[0]
        number_measured_obs = k[1]
        number_measured_all = "?"
        number_unique_obs = k[2]
        pdbx_rejected = str(int(k[1]) - int(k[7]))
        Rmerge_I_obs = p2f(k[5])
        meanI_over_sigI_obs = k[8]
        pdbx_Rsym_value = "?"
        pdbx_chi_squared = "?"
        pdbx_redundancy = "{0:.2f}".format(float(k[1]) / float(k[3]))
        percent_possible_obs = "?"
        pdbx_netI_over_sigmaI_obs = "?"
        number_possible = "?"
        number_unique_all = k[3]
        Rmerge_F_all = "?"
        Rmerge_F_obs = "?"
        Rmerge_I_all = "?"
        meanI_over_sigI_all = "?"
        percent_possible_all = p2f(k[4])
        pdbx_Rrim_I_all = p2f(k[9])
        pdbx_Rpim_I_all = "?"
        pdbx_number_anomalous = int(k[2]) - int(k[-1])
        pdbx_CC_half = p2f(k[-4].replace("*", "%"))

        table_line = (
            f"{pdbx_diffrn_id} {pdbx_ordinal} {d_res_high} {d_res_low} {number_measured_obs}"
            + f" {number_measured_all} {number_unique_obs} {pdbx_rejected} {Rmerge_I_obs}"
            + f" {meanI_over_sigI_obs} {pdbx_Rsym_value} {pdbx_chi_squared} {pdbx_redundancy}"
            + f" {percent_possible_obs} {pdbx_netI_over_sigmaI_obs} {number_possible} {number_unique_all}"
            + f" {Rmerge_F_all} {Rmerge_F_obs} {Rmerge_I_all} {meanI_over_sigI_all} {percent_possible_all}"
            + f" {pdbx_Rrim_I_all} {pdbx_Rpim_I_all} {pdbx_CC_half} {pdbx_number_anomalous}"
        )
        output.append(table_line)
    return "#\n" + n_line + "\n".join(output[::-1])


def parse_struct_keywords(protein_function, keywords):
    _struct_keywords = f"""# 
_struct_keywords.entry_id   UNNAMED 
_struct_keywords.text      'FragMAX, FragMAXapp, {keywords}' 
_struct_keywords.pdbx_keywords {protein_function} """
    return _struct_keywords


def parse_citation(citation_title, journal_abrev, journal_vol, first_page, last_page, year, PUBMED_id, DOI):
    _citation = f"""#
_citation.id                        primary 
_citation.title                     "{citation_title}"
_citation.journal_abbrev            "{journal_abrev}"
_citation.journal_volume            {journal_vol}
_citation.page_first                {first_page}
_citation.page_last                 {last_page}
_citation.year                      {year}
_citation.pdbx_database_id_PubMed   {PUBMED_id}
_citation.pdbx_database_id_DOI      {DOI}
"""

    return _citation


def parse_struct(mmcif, protein, library):
    model_id = mmcif.split("-")[-1].split("_")[0]
    if "Apo" in model_id:
        model_id = model_id.replace("Apo", "")
        state = "ground state model"
    else:
        state = f"changed state model for fragment {library_full_name}"
    _struct = f"""# 
_struct.entry_id   UNNAMED 
_struct.title      
'PanDDA analysis group deposition -- {protein_full_name} {state} {model_id} '
"""

    return _struct


def parse_exptl_crystal(matthews_coef, density_percent_sol):

    _exptl_crystal = f"""#
_exptl_crystal.id                           1
_exptl_crystal.density_Matthews             {matthews_coef}
_exptl_crystal.density_percent_sol          {density_percent_sol}
"""

    return _exptl_crystal


def parse_diffrn_radiation():
    _diffrn_radiation = f"""# 
_diffrn_radiation.diffrn_id                        1 
_diffrn_radiation.wavelength_id                    1 
_diffrn_radiation.pdbx_diffrn_protocol             'SINGLE WAVELENGTH' 
_diffrn_radiation.monochromator                    'Si(111)' 
_diffrn_radiation.pdbx_monochromatic_or_laue_m_l   M 
"""
    return _diffrn_radiation


def parse_entity_src_nat(nat_strain, organism_scientific, taxon_ID):
    _entity_src_nat = f"""# 
_entity_src_nat.entity_id                  1 
_entity_src_nat.pdbx_src_id                1 
_entity_src_nat.pdbx_organism_scientific   {organism_scientific} 
_entity_src_nat.strain                     {nat_strain}
_entity_src_nat.pdbx_ncbi_taxonomy_id      {taxon_ID}
_entity_src_nat.details                    ? 
"""
    return _entity_src_nat


def parse_entity_src_gen(
    beg_seq_num,
    end_seq_num,
    gene_src_gene,
    organism_scientific,
    taxon_ID,
    host_org_scientific_name,
    host_org_ncbi_taxonomy_id,
):

    _entity_src_gen = f"""# 
_entity_src_gen.entity_id                          1 
_entity_src_gen.pdbx_src_id                        1 
_entity_src_gen.pdbx_alt_source_flag               sample 
_entity_src_gen.pdbx_seq_type                      'Biological sequence' 
_entity_src_gen.pdbx_beg_seq_num                   {beg_seq_num} 
_entity_src_gen.pdbx_end_seq_num                   {end_seq_num} 
_entity_src_gen.gene_src_common_name               ? 
_entity_src_gen.gene_src_genus                     ? 
_entity_src_gen.pdbx_gene_src_gene                 {gene_src_gene} 
_entity_src_gen.gene_src_species                   ? 
_entity_src_gen.gene_src_strain                    ? 
_entity_src_gen.gene_src_tissue                    ? 
_entity_src_gen.gene_src_tissue_fraction           ? 
_entity_src_gen.gene_src_details                   ? 
_entity_src_gen.pdbx_gene_src_fragment             ? 
_entity_src_gen.pdbx_gene_src_scientific_name      {organism_scientific} 
_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id     {taxon_ID} 
_entity_src_gen.pdbx_gene_src_variant              ? 
_entity_src_gen.pdbx_gene_src_cell_line            ? 
_entity_src_gen.pdbx_gene_src_atcc                 ? 
_entity_src_gen.pdbx_gene_src_organ                ? 
_entity_src_gen.pdbx_gene_src_organelle            ? 
_entity_src_gen.pdbx_gene_src_cell                 ? 
_entity_src_gen.pdbx_gene_src_cellular_location    ? 
_entity_src_gen.host_org_common_name               ? 
_entity_src_gen.pdbx_host_org_scientific_name      {host_org_scientific_name} 
_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id     {host_org_ncbi_taxonomy_id} 
_entity_src_gen.host_org_genus                     ? 
_entity_src_gen.pdbx_host_org_gene                 ? 
_entity_src_gen.pdbx_host_org_organ                ? 
_entity_src_gen.host_org_species                   ? 
_entity_src_gen.pdbx_host_org_tissue               ? 
_entity_src_gen.pdbx_host_org_tissue_fraction      ? 
_entity_src_gen.pdbx_host_org_strain               ? 
_entity_src_gen.pdbx_host_org_variant              ? 
_entity_src_gen.pdbx_host_org_cell_line            ? 
_entity_src_gen.pdbx_host_org_atcc                 ? 
_entity_src_gen.pdbx_host_org_culture_collection   ? 
_entity_src_gen.pdbx_host_org_cell                 ? 
_entity_src_gen.pdbx_host_org_organelle            ? 
_entity_src_gen.pdbx_host_org_cellular_location    ? 
_entity_src_gen.pdbx_host_org_vector_type          ? 
_entity_src_gen.pdbx_host_org_vector               ? 
_entity_src_gen.host_org_details                   ? 
_entity_src_gen.expression_system_id               ? 
_entity_src_gen.plasmid_name                       ? 
_entity_src_gen.plasmid_details                    ? 
_entity_src_gen.pdbx_description                   ? 
"""
    return _entity_src_gen


def p2f(x):
    return "{0:.3f}".format(float(x.strip("%")) / 100)


def new_line_list(line):
    # when processing logfiles, rearrange
    return [x.replace(":", "").replace("|", "") for x in line.split()[1:]]


def write_to_mmcif(mmcif, content):
    if type(content) is list:
        print("content parsed from list")
        content = "".join(content)
    else:
        print("content parsed from string")
    with open(mmcif, "w") as w:
        w.write(content)
        print(mmcif + " wrote succesfully")


def get_value_field(mmcif, field):
    with open(mmcif, "r") as r:
        c = r.readlines()
    for line in c:
        if field + " " in line:
            return line


def replace_line(mmcif, field, new_value):
    new_value = str(new_value)
    with open(mmcif, "r") as r:
        c = r.readlines()
    if field + " " in "".join(c):
        for n, line in enumerate(c):

            if field + " " in line:
                print("Field found in " + mmcif)
                old_value = line.split()[-1]

                if new_value != old_value:
                    print("Changing: " + old_value + " to " + new_value)
                    print("for " + field)
                    c[n] = c[n].replace(old_value, new_value)

                else:
                    print("New value is equal to old value for " + mmcif)
                    print("# Skipping")
        with open(mmcif, "w") as w:
            w.write("".join(c))
    else:
        print("Entry " + field + " not found in " + mmcif)
        print("# Skipping")


def modify_section_mmcif(mmcif, section, mode, new_section=""):
    print("Making changes to " + mmcif)
    with open(mmcif, "r") as r:
        c = r.readlines()
    eventLines = [n for n, line in enumerate(c) if section + "." in line]
    if len(eventLines) != 0:
        print("Section " + section + " found!")
        bangs = [n for n, line in enumerate(c) if "#" in line]

        first = min(bangs, key=lambda x: abs(x - eventLines[0]))
        if first == bangs[-1]:
            last = first
        else:
            last = bangs[bangs.index(first) + 1]
        print(first, last)
    else:
        print("Section " + section + " not in " + mmcif)
    if mode == "delete":
        cropFile = c[:first] + c[last + 1 :]
        removed = c[first:last]
        with open(mmcif, "w") as w:
            w.write("".join(cropFile))
        # return cropFile, removed
    elif mode == "replace":
        if type(new_section) is not list:
            new_section = [x + "\n" for x in new_section.splitlines()]
        print("replacing sections")
        print("--- OLD SECTION ---")
        print("".join(c[first:last]))
        print("--- NEW SECTION ---")
        print("".join(new_section))
        newFile = c[:first] + new_section + ["#\n"] + c[last + 1 :]
        with open(mmcif, "w") as w:
            w.write("".join(newFile))
        # return newFile, new_section
    elif mode == "check":
        print("Checking for a section")
        if section in "".join(c):
            print("--- OLD SECTION ---")
            print("".join(c[first:last]))
        else:
            print("Section " + section + " not found in " + mmcif)
    elif mode == "add":
        print("Adding new section")
        print(section)
        if section in "".join(c):
            print("Section is already present. Choose mode==replace to replace it.")
            print("Skipping")

        else:
            if type(new_section) is not list:
                new_section = [x + "\n" for x in new_section.splitlines()]
            if not "\n" in c[-1]:
                c.append("\n")
            newFile = c + new_section + ["#\n"]
            with open(mmcif, "w") as w:
                w.write("".join(newFile))
    else:
        print("Section not found in " + mmcif)


def select_search_file(dataset_file):
    with open(dataset_file, "r") as readFile:
        pdb_content = readFile.readlines()
    use_pdb = False
    for line in pdb_content:
        if "FIT TO DATA USED IN REFINEMENT (IN BINS)." in line:
            use_pdb = True
    return use_pdb


def find_refine_table_phenix(file_content, file_format, logFile):
    if file_format == "pdb":
        for n, line in enumerate(file_content):
            if "FIT TO DATA USED IN REFINEMENT (IN BINS)." in line:
                table_ref = file_content[n : n + 20]
                # end_of_table=table_ref.index('REMARK   3  \n')
                for n, x in enumerate(table_ref):
                    if len(x.split()) == 2:
                        end_of_table = n
                table_ref = [x.split()[2:] for x in table_ref[2:end_of_table]][::-1]
                return table_ref
    if file_format == "reflog":
        for n, line in enumerate(file_content):
            if "| Bin     Resolution   Compl.  No. Refl.    R-factors          Targets        |" in line:
                table = file_content[n + 2 : n + 20]
                for n, i in enumerate(table):
                    if "|--------------------------" in i:
                        k = n
                        table = table[:k]
        return [new_line_list(x) for x in table][::-1]


def find_refine_table_refmac(file_content, file_format, logFile):
    for line in file_content:
        binSize = "1"
        if "Resolution limits" in line:
            low, high = line.split()[-2:]
        if "Percentage observed" in line:
            completeness = float(line.split()[-1]) / 100
        if "Number of used reflections" in line:
            number_reflns_all = line.split()[-1]
        if "Percentage of free reflections" in line:
            free_perc = line.split()[-1]
            if float(free_perc) > 1.0:
                free_perc = float(free_perc) / 100
            elif float(free_perc) == 0.0:
                free_perc = 0.0

            number_reflns_R_work = int((1.0 - free_perc) * float(number_reflns_all))
            number_reflns_R_free = int(number_reflns_all) - number_reflns_R_work
    with open(logFile.replace(".log", ".pdb"), "r") as readFile:
        ref_pdb = readFile.readlines()
    for line in ref_pdb:
        if "R VALUE     (WORKING + TEST SET)" in line:
            ls_R_factor_R_free = line.split()[-1]
        if "R VALUE            (WORKING SET)" in line:
            ls_R_factor_R_work = line.split()[-1]
    t = [
        "1",
        low,
        "-",
        high,
        completeness,
        number_reflns_R_work,
        number_reflns_R_free,
        ls_R_factor_R_work,
        ls_R_factor_R_free,
    ]
    return [[str(x) for x in t]]


def fix_entry(mmcif, path):
    dataset = mmcif.split("/")[-1].replace(".mmcif", "")

    mmcif = f"{path}/{dataset}.mmcif"
    pdbFile = f"{path}/{dataset}.pdb"
    logFile = f"{path}/{dataset}.log"
    correctFile = f"{path}/{dataset}_CORRECT.LP"
    unmrgFile = f"{path}/{dataset}_data_F.mtz"

    write_to_mmcif(mmcif, fix_occupancies_ensemmble(mmcif))

    write_to_mmcif(mmcif, fix_occupancies_conformer(mmcif))

    modify_section_mmcif(mmcif, "_reflns_shell", "add", retrive_reflns_shell_section(correctFile))
    modify_section_mmcif(mmcif, "_pdbx_deposit_group", "add", create_group_description(protein, library, group_dep_id))
    modify_section_mmcif(
        mmcif,
        "_entity_src_gen",
        "add",
        parse_entity_src_gen(
            beg_seq_num,
            end_seq_num,
            gene_src_gene,
            organism_scientific,
            taxon_ID,
            host_org_scientific_name,
            host_org_ncbi_taxonomy_id,
        ),
    )

    modify_section_mmcif(
        mmcif, "_refine_ls_shell", "replace", new_section=generate_refine_ls_shell_table(pdbFile, logFile)
    )
    modify_section_mmcif(mmcif, "_audit_author", "replace", new_section=parse_audit_author(u))
    modify_section_mmcif(mmcif, "_pdbx_contact_author", "replace", new_section=create_pdbx_contact_author(u))
    modify_section_mmcif(mmcif, "_citation_author", "replace", new_section=parse_citation_author(u))
    modify_section_mmcif(mmcif, "_pdbx_entity_instance_feature", "add", parse_pdbx_entity_instance_feature())
    modify_section_mmcif(mmcif, "_software", "replace", _software)
    modify_section_mmcif(mmcif, "_struct_keywords", "replace", parse_struct_keywords(protein_function, keywords))
    modify_section_mmcif(
        mmcif,
        "_citation",
        "replace",
        parse_citation(citation_title, journal_abrev, journal_vol, first_page, last_page, year, PUBMED_id, DOI),
    )
    modify_section_mmcif(mmcif, "_struct", "replace", parse_struct(mmcif, protein, library))
    modify_section_mmcif(mmcif, "_diffrn_radiation", "replace", parse_diffrn_radiation())
    modify_section_mmcif(
        mmcif, "_entity_src_nat", "replace", parse_entity_src_nat(nat_strain, organism_scientific, taxon_ID)
    )

    replace_line(mmcif, "_refine.pdbx_method_to_determine_struct", "'FOURIER SYNTHESIS'")
    replace_line(mmcif, "_diffrn_detector.pdbx_collection_date", collection_date)

    for k, v in parse_statistics_from_logs(logFile, correctFile).items():
        replace_line(mmcif, k, v)


def generate_index(export_datasets):
    # Create index file for group dep based on exported mmCIF files

    def new_entry(model):
        entry = f"""label: {model}\n""" f"""model: {model}.mmcif\n""" f"""sf: {model}_sf.mmcif\n"""
        return entry

    pdbFile = export_datasets[0]
    path = "/".join(pdbFile.split("/")[:-1])
    protein = pdbFile.split("/")[-1].split("-")[0]
    library = pdbFile.split("/")[-1].split("-")[1]

    with open(f"{path}/index.txt", "w") as writeFile:
        for dt in export_datasets:
            pdbFile = dt
            run = dt.split("-")[-1].split("_")[-1].split(".")[0]
            if "Apo" in dt:
                fragment = dt.split("-")[-1].split("_")[0]
                sub = ""
            else:
                fragment = dt.split("-")[-1].split("_")[0][:-1]
                sub = dt.split("-")[-1].split("_")[0][-1]
            model = f"{protein}-{library}-{fragment}{sub}_{run}"
            print(pdbFile, fragment, sub, run)
            if os.path.exists(f"{path}/{model}.mtz"):
                writeFile.write(new_entry(model))
                generate_log_script(pdbFile, fragment, sub, run)
                generate_data_template(pdbFile, fragment, sub, run)
                print("done")
                cmd = f"pdb_extract -ext log_script_{fragment}{sub}_{run}.inp "
                print(cmd)
                subprocess.call(cmd, shell=True)
                mmcif = f"{path}/{model}.mmcif"
                print(mmcif)
                fix_entry(mmcif, path)


dataset_export_list = glob(f"{path}/*pdb")
print(str(len(dataset_export_list)) + " datasets found.")
print([x for x in dataset_export_list])
generate_index(dataset_export_list)


if os.path.exists(path + "/SF_xxxx__.mmcif"):
    os.remove(path + "/SF_xxxx__.mmcif")


subprocess.call("zip group_python.zip *mmcif index.txt", shell=True)
