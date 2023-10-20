import os
import time

from constants.constants import Televir_Metadata_Constants

entrez_direct_bin = os.path.join(
    Televir_Metadata_Constants.BINARIES["ROOT"],
    Televir_Metadata_Constants.BINARIES["software"]["entrez_direct"],
    "bin",
)


def query_taxid(description: str, tmp_dir) -> list:
    # 20 second max wait time
    wait_time_max = 20
    print(description)

    tempfilename = (
        "query_"
        + description.replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace(";", "")
        + ".txt"
    )
    tempfilename = os.path.join(tmp_dir, tempfilename)
    tmp_bash_script = os.path.join(tmp_dir, "tempfile.sh")

    query = f'{entrez_direct_bin}/esearch -db genome -query "{description}" | {entrez_direct_bin}/efetch -format docsum | {entrez_direct_bin}/xtract -pattern DocumentSummary -element TaxId > {tempfilename}'
    with open(tmp_bash_script, "w") as tempfile:
        tempfile.write(query)

    # deploy job in background
    open(tempfilename, "w").close()

    os.system(f"chmod +x {tmp_bash_script}")
    os.system(f"sh {tmp_bash_script} &")
    start_time = time.time()
    # while file is empty and time is less than 20 seconds:

    while not os.path.getsize(tempfilename) > 0:
        time.sleep(0.1)
        if time.time() - start_time > wait_time_max:
            break

    if os.path.exists(tempfilename):
        tempfile = open(tempfilename, "r")
        stdout = tempfile.read()
        tempfile.close()
        os.system(f"rm {tempfilename}")

        stdout = stdout.strip().split("\n")
        if stdout[0] == "":
            return ["NA"]
        return stdout
    else:
        os.remove(tempfilename)
        print("NA")
        return []


def entrez_fetch_taxid_from_org_description(
    description: str, db: str = "taxonomy"
) -> list:
    from Bio import Entrez

    Entrez.email = "joao.dourado@insa.min-saude.pt"
    try:
        handle = Entrez.esearch(db="genome", term=description)
        record = Entrez.read(handle)
        handle.close()
    except:
        return ["NA"]
    # print(record)
    # print(record)
    # print(description)
    print(description)
    if len(record["IdList"]) == 0:
        return ["NA"]
    return record["IdList"]


def remove_parenthesis(description: str) -> str:
    if "(" in description:
        return description.split("(")[0].strip()
    return description


def cut_after_virus(description: str) -> str:
    if "virus" in description:
        desc = description.split("virus")[0] + "virus"
        # check if space before "virus"
        return desc
    return description


def remove_last_word(description: str) -> str:
    ## keep first two words
    return " ".join(description.split(" ")[:2]).strip()


def taxid_passes_test(taxid_list: list) -> bool:
    if taxid_list[0] == "NA":
        return False

    return True


def specific_amends(description: str) -> str:
    if "Gripe" in description:
        return "Influenza A"
    if "Rinovírus" in description:
        return description.replace("Rinovírus", "Rhinovirus")

    if "Human herpesvirus 6 (HHV-6)" in description:
        return "Human herpesvirus 6"

    return description


def entrez_fetch_taxid_from_org_description_curate(
    description: str, tmp_dir: str = "/tmp"
) -> list:
    description = specific_amends(description)

    taxid = query_taxid(description, tmp_dir)

    if not taxid_passes_test(taxid):
        desc_small = remove_parenthesis(description)
        desc_small = cut_after_virus(desc_small)
        taxid = query_taxid(desc_small, tmp_dir)

    if not taxid_passes_test(taxid):
        taxid = entrez_fetch_taxid_from_org_description(description)

    if not taxid_passes_test(taxid):
        desc_small = remove_parenthesis(description)
        desc_small = cut_after_virus(desc_small)
        taxid = entrez_fetch_taxid_from_org_description(desc_small)

    if not taxid_passes_test(taxid) and len(description.split(" ")) > 2:
        desc_small = remove_last_word(description)
        taxid = entrez_fetch_taxid_from_org_description(desc_small)

    return taxid
