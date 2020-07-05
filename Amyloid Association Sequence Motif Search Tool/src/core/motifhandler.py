from Bio import ExPASy
from Bio import SeqIO
import re
import tkinter as tk
import tkinter.ttk as ttk
import requests as r
from io import StringIO


class MotifHandler:
    def getSequence(self, code):
        uniMotif = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
        if re.fullmatch(uniMotif, code):
            uniURL = "http://www.uniprot.org/uniprot/"
            fastaURL = uniURL + code + ".fasta"
            response = r.post(fastaURL)
            fastaText = ''.join(response.text)
            seq = StringIO(fastaText)
            parsedSeq = list(SeqIO.parse(seq, 'fasta'))
            return parsedSeq
        else:
            return None
    
    def queryMotif(self, motif, sequence):
        search = re.finditer(motif, str(sequence))
        return list(search)
