from Bio import Entrez

def count_occurrences(tool_name):
    Entrez.email = "your.email@example.com"  # Please enter your email address here
    handle = Entrez.esearch(db="pubmed", term=tool_name)
    record = Entrez.read(handle)
    return int(record["Count"])

## pipe list of tools, or topics to count citations
tools = [sys.argv[1]]

for tool in tools:
    count = count_occurrences(tool)
    print(f"The tool {tool} is mentioned in {count} articles in PubMed.")
