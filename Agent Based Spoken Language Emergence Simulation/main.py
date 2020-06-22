from Agent import Agent
import random
import statistics

##VARIABLES##
agentNum = 2
lsRatio = 1
streakbreak = 25
attributes = ["big","small","round","square","long","short","red","yellow",
                  "green","dark","light","sharp","blunt","hot","cold","clear",
                  "opaque","metallic","hairy","heavy","light","dense","hollow",
                  "brown","blue","black","white"]
symbols = ['æ', 'ɒ', 'ɛ', 'ɪ', 'ʊ', 'ʌ', 'ə', 'i', 'u', 'ɑ', 'b',
          'd', 'ð', 'f', 'ɡ', 'h', 'j', 'k', 'l', 'm', 'n', 'ŋ',
          'p', 'r', 's', 'ʃ', 't', 'θ', 'v', 'w', 'z', 'ʒ']
agents = [None]*agentNum
rounds = 1000000
speakers = []
listeners = []
streak = 0
i = 0

seed = "seed"
random.seed(seed)
##---------##

def createAgents():
    global agents
    for i in range(agentNum):
        agent = Agent(attributes,symbols)
        agents[i] = agent

def setAgents():
    global speakers
    global listeners
    speakers = random.sample(agents, int(len(agents)/(lsRatio+1)))
    listeners = [x for x in agents if x not in speakers]

def setRoles():
    step = int(len(listeners)/len(speakers))
    rem  = int(len(listeners)%len(speakers))
    count = 0
    for speaker in speakers:
        listenList = []
        for i in range(step):
            listenList.append(listeners[count+i])
        count = count + step
        speaker.listeners = listenList
    for i in range(rem):
        speakers[i].listeners.append(listeners[count])
        count = count + 1

def assignJobs():
    for speaker in speakers:
        for listener in speaker.listeners:
            i = random.randint(1,7)
            listener.job = random.sample(attributes,i)


            
def main():
    print("Start: ", seed)
    print("Population size: ", agentNum)
    print("Listener/speaker ratio: ", lsRatio)
    print("Streak Limit: ", streakbreak)
    print("Attributes: ", len(attributes))
    createAgents()
    for i in range(rounds):
        setAgents()
        setRoles()
        assignJobs()
        correct = 0
        incorrect = 0
        for speaker in speakers:
            for listener in speaker.listeners:
                result = listener.listen(speaker.speak(listener.job))
                if set(listener.job) == set(result):
                    correct += 1
                else:
                    incorrect += 1
                listener.lupdate(result)
                speaker.supdate(listener.job, result)
        if incorrect == 0:
            streak += 1
            if streak >= streakbreak:
                break
        else:
            streak = 0
    print("Done Round: ", i)
    print(agents[0].memory)

if __name__ == "__main__":
    main()
