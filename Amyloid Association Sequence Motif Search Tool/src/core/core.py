from Bio import ExPASy
from Bio import SeqIO
import re
import tkinter as tk
import tkinter.ttk as ttk
import requests as r
from io import StringIO

from motifhandler import MotifHandler 
from infoviewer import InfoViewer

# Base tkinter app. Allows for searching of UNIPROT accession code's sequence against the motif
# Displays matches of the motif for each code and allows for user to choose any proteins to
# display more information about.
class MainApp(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.master = master
        self.setConfigs()
        self.createWidgets()
        # Prosite motif into regex form, with a lookbehind to match overlapping match
        self.motif = r'(?=([^P][^PKRHW][VLSWFNQ][ILTYWFN][FIY][^PKRH]))'
        # Handles the retrieval of sequence data from UNIPROT codes.
        self.motifHandler = MotifHandler()
        self.matches = []
        self.sequences = []
        self.codes = []

    # Set configuration settings 
    def setConfigs(self):
        self.master.title("Motif Searcher")
        self.master.resizable(False, False)

    # Set up the window's UI
    def createWidgets(self):
        #
        self.titleLabel = tk.Label(self, text="Searching for motif: {P}-{PKRHW}-[VLSWFNQ]-[ILTYWFN]-[FIY]-{PKRH}")
        self.titleLabel.grid(row=0,column=0,columnspan=5)
        #
        self.codeLabel = tk.Label(self, text="Code(s):")
        self.codeLabel.grid(row=1,column=0, sticky="n")
        #
        self.codeText = tk.Text(self, height=8, width=30)
        self.codeText.grid(row=1, column=1, sticky="ew")
        self.codeText.grid_columnconfigure(3, weight=1)
        codeTextScrollY = tk.Scrollbar(self, command=self.codeText.yview)
        codeTextScrollY.grid(row=1, column=2, stick='nsw')
        self.codeText['yscrollcommand'] = codeTextScrollY.set
        self.codeTextOut(None)
        #
        self.enterCodeButton = tk.Button(self, text="Enter", command=self.enterButton)
        self.enterCodeButton.grid(row=1, column=3, columnspan=1, sticky="w")
        #
        self.displayText = tk.Text(self, width=80, height=15, pady=5, state="disabled")
        self.displayText.grid(row=2, column=0, columnspan=4, pady=(5,0), sticky="ew")
        displayTextScrollY = tk.Scrollbar(self, command=self.displayText.yview)
        displayTextScrollY.grid(row=2, column=4, stick='nse', pady=(5,5))
        self.displayText['yscrollcommand'] = displayTextScrollY.set
        #
        self.codeText.bind("<FocusIn>", self.codeTextIn)
        self.codeText.bind("<FocusOut>", self.codeTextOut)
        #
        self.viewMoreLabel = tk.Label(self, text="Enter the line numbers of the motifs you would like to view more information about below:",)
        self.viewMoreLabel.grid(row=3, column=1)
        #
        self.viewMoreFrame = tk.Frame(self)
        self.viewMoreFrame.grid(row=4, column=0, columnspan=5)
        self.viewMoreEnter = tk.Entry(self.viewMoreFrame)
        self.viewMoreEnter.grid(row=0, column=0, sticky="ew")
        #
        self.viewMoreButton = tk.Button(self.viewMoreFrame, text="Enter", command=self.moreButton)
        self.viewMoreButton.grid(row=0, column=1, padx=(5,0))

    # Handle logic of user changing focus to the code entry widget
    # If the example text was entered the box will clear when entered
    def codeTextIn(self, event):
        hintText = "Eg:\nP10636\nP10997\nP04156\nP05067\nP02647"
        if (self.codeText.get(1.0,"end-1c") == hintText):
            self.codeText.delete(1.0,"end")

    # Handle logic of user changing focus away from the code entry widget
    # If there is no text in the widget the example text will re-fill the widget
    def codeTextOut(self, event):
        hintText = "Eg:\nP10636\nP10997\nP04156\nP05067\nP02647"
        if (self.codeText.get(1.0,"end-1c") == ""):
            self.codeText.insert(tk.INSERT, hintText)

    # Loop through codeText, sep bu newline to get each code
    # If the line isnt a valid code give a message in the output box
    # Handle actual codes accordingly
    def enterButton(self):
        self.codes = self.codeText.get(1.0,"end").splitlines()
        self.displayText.config(state="normal")
        self.displayText.delete(1.0,"end")
        self.displayText.config(state="disabled")
        lineNum = 0
        self.sequences = []
        self.matches = []
        for code in self.codes:
            lineNum += 1
            self.displayText.config(state="normal")
            self.displayText.insert(tk.INSERT, "" + str(self.handleCode(code.upper(), lineNum)) + "\n")
        self.displayText.config(state="disabled")

    # Create the line of text shown in the display for a given Uniprot code
    # Also appends sequences and matches to an attribute list
    def handleCode(self, code, num):
        holding = self.motifHandler.getSequence(code)
        if holding:
            self.sequences.append(holding)
            for sequence in self.sequences[num-1]:
                output = str(num) + ". " + code + ":  "
                matches = self.motifHandler.queryMotif(self.motif, sequence.seq)
                self.matches.append(matches)
                for match in matches:
                    output = output + "[" + str(match.start(1)+1) + ", " + str(match.end(1)) + "]:" + match.group(1) + " "
                return output
        else:
            self.sequences.append(None)
            self.matches.append(None)
            return str(num) + ". Not a valid UniProt accession number. Please check your format and try again." 

    # Check the number entered into the textbox is a valid line number
    # If not create an error message telling the user to correct it
    # If it is create a window displaying more information about the match
    def moreButton(self):
        texts = self.viewMoreEnter.get().split(" ")
        nums = []
        # Turn viewMore entry into line number integers
        for text in texts:
            try:
                nums.append(int(text))
            except:
                nums.append(None)
        for num in nums:
            try:
                if self.sequences[num-1] and self.codes[num-1] and self.matches[num-1]:
                    viewmore = InfoViewer(tk.Tk(), self.codes[num-1], self.sequences[num-1], self.matches[num-1], num)
                    viewmore.pack()
                elif(not self.codes[num-1]):
                    print("No codes available for line " + str(num))
                elif(not self.sequences[num-1]):
                    print("No sequences available for line " + str(num))
                elif(not self.matches[num-1]):
                    print("No matches available for line " + str(num))
                else:
                    print("No information can be displayed for line "+ str(num))
            except Exception as e:
                print("Line number entered is errenuous. Try again.")
                print(e)

def main():
    window = MainApp(tk.Tk())
    window.config(pady=5, padx=5)
    window.pack(expand = True)
    window.mainloop()

