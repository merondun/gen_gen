import os
import sys
import ftplib
from Bio import Entrez

def download_refseq_files(refseq_accession):
    Entrez.email = "your_email@example.com"  # Replace with your email address

    # Retrieve assembly summary
    handle = Entrez.esearch(db="assembly", term=refseq_accession)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        print(f"No assembly found for RefSeq accession: {refseq_accession}")
        return

    assembly_id = record["IdList"][0]

    # Fetch assembly report
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    summary = Entrez.read(handle)
    handle.close()

    ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]

    if not ftp_path:
        print(f"No FTP path found for RefSeq accession: {refseq_accession}")
        return

    ftp_server = "ftp.ncbi.nlm.nih.gov"
    ftp_dir = "/".join(ftp_path.split("/")[3:])

    # Download files
    dir_name = refseq_accession
    os.makedirs(dir_name, exist_ok=True)

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
    if len(sys.argv) != 2:
        print("Usage: python download_refseq_files.py <RefSeq_accession>")
        sys.exit(1)

    refseq_accession = sys.argv[1]
    download_refseq_files(refseq_accession)
