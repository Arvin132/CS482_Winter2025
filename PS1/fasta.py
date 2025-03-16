

def read_fasta_file(file: str) -> list[tuple[str, str]]:
    
    """
    Reads a FASTA file and extracts the headers and sequences.

    Parameters:
    -----------
    file : str
        Path to the FASTA file to be read.

    Returns:
    --------
    list[tuple[str, str]]
        Each tuple contains:
        - header
        - The corresponding sequence
    
    Raises:
    -------
    ValueError
        If the file is not in a valid FASTA format (e.g., missing '>' at the start).
    FileNotFoundError
        If the specified file does not exist.

    Notes:
    ------
    - Assumes the FASTA file does not contain any multi-line headers.
    - The sequence is treated as a single string without newlines.
    """
    
    retVal = []
    
    with open(file, mode='r+') as f:
        ch = ch = f.read(1)
        while(True): 
            if(ch == ""): break
            if (ch != '>'): raise ValueError("Bad Input FASTA file")
            
            header = ""
            ch = f.read(1)
            while(ch != "\n"):
                header += ch
                ch = f.read(1)
                
            content = ""
            ch = f.read(1)
            while(ch != ">" and ch != ""):
                if (ch != '\n'):
                    content += ch
                ch = f.read(1)
            
            retVal.append((header, content))
    
    return retVal
