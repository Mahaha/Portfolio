import collections
import random
#import Levenshtein as le

class Agent:


    def __init__(self, atts, symbs):
        self.memory = {}
        self.listeners = []
        self.assos = {}
        self.atts = atts
        self.symbols = symbs
        self.delta = 0.1

    def addMem(self, signal, attribute, score):
        try:
            self.memory[signal][attribute] = score
        except:
            self.memory[signal] = {}
            self.memory[signal][attribute] = score
                    
    def cleanup(self):
        sigs = []
        attrs = {}
        for signal, atts in self.memory.items():
            for attr, scr in atts.items():
                if scr <= 0:
                    attrs[signal] = attr
                if scr > 1:
                    self.memory[signal][attr] = 1
        for sig, att in attrs.items():
            self.memory[sig].pop(att)
        for signal, atts in self.memory.items():
            if len(atts) <= 0:
                sigs.append(signal)
        for sig in sigs:
            self.memory.pop(sig)
            
    def speak(self, attributes):
        signal = ""
        for attribute in attributes:
            word = ""
            score = 0
            for sig, atts in self.memory.items():
                for att, scr in atts.items():
                  if att == attribute and scr > score:
                      word = sig
                      score = scr
            if word == "":
                k = random.randint(1,8)
                word = word.join(random.sample(self.symbols, k))
                self.addMem(word, attribute, self.delta)
            self.assos[word] = attribute
            signal = signal + word + " "
        signal = signal.rstrip()
        return signal

    def listen(self, signalin):
        signals = signalin.split(" ")
        attributes = []
        for signal in signals:
            closest = None
            dist = None
            for sig in self.memory.keys():
                if signal == sig:
                    closest = sig
                    dist = 0
            highest = random.choice(self.atts)
            score = 0
            if closest != None:
                for att, scr in self.memory[closest].items():
                    if highest == None or scr > score:
                        highest = att
                        score = scr
            self.assos[signal] = highest
            attributes.append(highest)
        return attributes
      
    def supdate(self, job, result):
        corrs = [i for i in result if i in job]
        incorrs = [i for i in result if i not in job]
        for sig, att in self.assos.items():
            if att in corrs:
                self.memory[sig][att] += self.delta
                for signal, attrs in self.memory.items():
                    for attribute, scr in attrs.items():
                        if signal != sig and attribute == att:
                            self.memory[signal][attribute] -= self.delta
            else:
                self.memory[sig][att] -= self.delta
        self.assos = {}
        self.cleanup()

    def lupdate(self, result):
        job = self.job
        corrs = [i for i in result if i in job]
        incorrs = [i for i in result if i not in job]
        
        for sig, att in self.assos.items():
            
            if att in corrs:
                try:
                    self.memory[sig][att] += self.delta*2
                except:
                    self.addMem(sig,att,self.delta*2)
                for signal, attrs in self.memory.items():
                    for attribute, scr in attrs.items():
                        if signal == sig and attribute != att:
                            self.memory[signal][attribute] -= self.delta
            else:
                try:
                    self.memory[sig][att] -= self.delta
                except:
                    pass
        self.assos = {}
        self.cleanup()
