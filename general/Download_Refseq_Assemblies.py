import os
import sys
import ftplib
from Bio import Entrez

# Submit with refseq ID as a positional argument 
def download_refseq_files(refseq_accession):
    # Set your email to communicate with NCBI
    Entrez.email = "heritabiltiies@gmail.com"  #please use your own   

    # Perform a search on the NCBI Assembly database using the provided refseq_accession
    handle = Entrez.esearch(db="assembly", term=refseq_accession)
    record = Entrez.read(handle)
    handle.close()

    # If no records are found for the provided refseq_accession, the function will terminate here.
    if not record["IdList"]:
        print(f"No assembly found for RefSeq accession: {refseq_accession}")
        return

    assembly_id = record["IdList"][0]

    # Retrieve the summary of the assembly with the identified id
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    summary = Entrez.read(handle)
    handle.close()

    # Extract the FTP path from the summary 
    ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]

    # If the ftp path is not found, the function will terminate here
    if not ftp_path:
        print(f"No FTP path found for RefSeq accession: {refseq_accession}")
        return

    ftp_server = "ftp.ncbi.nlm.nih.gov"
    ftp_dir = "/".join(ftp_path.split("/")[3:])

    # Create a new directory with the name of the refseq_accession if it does not exist
    dir_name = refseq_accession
    os.makedirs(dir_name, exist_ok=True)

    # Establish a connection to the ftp server and download the files with the extensions .gz and .txt
    try:
        ftp = ftplib.FTP(ftp_server)
        ftp.login()
        ftp.cwd(ftp_dir)

        for file_name in ftp.nlst():
            if file_name.endswith((".gz", ".txt")):
                print(f"Downloading file: {file_name}")
                local_file = os.path.join(dir_name, file_name)
                with open(local_file, "wb") as f:
                    ftp.retrbinary(f"RETR {file_name}", f.write)

        ftp.quit()
    except Exception as e:
        print(f"Error downloading files: {e}")


if __name__ == "__main__":
    # Check if the command line argument (RefSeq accession) is provided.
    if len(sys.argv) != 2:
        print("Usage: python download_refseq_files.py <RefSeq_accession>")
        sys.exit(1)

    # Call the function with the command line argument
    refseq_accession = sys.argv[1]
    download_refseq_files(refseq_accession)


