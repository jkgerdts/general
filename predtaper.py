def predtaper(doses,weeks):
    #def dose_sig(dos):
        
    dose_split = {80:[20,20,20,20],60:[20,20,20],50:[20,20,10],40:[20,20],35:[20,10,5],30:[20,10],25:[20,5],20:[20],15:[10,5],12.5:[10,2.5],10:[10],7.5:[5,2.5],5:[5],2.5:[2.5]}
    sig = {80:'80mg (20x4)',60:'60mg (20x3)',50:'50mg (20+20+10)',40:'40mg (20+20)',35:'35mg (20+10+5)',30:'30mg (20+10)',25:'25mg (20+5)',20:'20mg',15:'15mg (10+5)',12.5:'12.5mg (10+2.5)',10:'10mg',7.5:'7.5mg (5+2.5)',5:'5mg',2.5:'2.5mg'}
    import datetime as dt
    dates = [dt.datetime.now()+dt.timedelta(days=1)]
    running = dates[0]
    runningweek = 0
    intervals = []
    daycounts = []
    weekranges = []
    def weekrangestring(weekrange):
        if weekrange[1]-weekrange[0]>0:
            return 'Weeks {}-{}: '.format(weekrange[0],weekrange[1])
        else:
            return 'Week {}: '.format(weekrange[0])
    for w in weeks:
        weekranges.append([runningweek+1,runningweek+w])
        runningweek = runningweek+w
        stop = running+dt.timedelta(days=7*w)
        intervals.append({'start':running,'stop':stop})
        interval = stop-running
        daycounts.append(interval.days)
        running = stop+dt.timedelta(days=1)
    pills = []
    
    # loop: for each dose, for each sub-dose, for each day, add the sub-dose to running list
    for i in range(len(doses)):
        for d in dose_split[doses[i]]:
            for j in range(daycounts[i]):
                pills.append(d)
                
    pillsU = list(set(pills))
    pillsC = []
    for pu in pillsU:
        pillsC.append(len([x for x in pills if x==pu]))
    
    for I,dos,dayct,wr in zip(intervals,doses,daycounts,weekranges):
        print(weekrangestring(wr)+'('+I['start'].strftime("%m/%d/%y")+'-'+I['stop'].strftime("%m/%d/%y")+'): ',sig[dos])
    for dos,wr in zip(doses,weekranges):
        print(weekrangestring(wr),sig[dos])
    for pu,pc in zip(pillsU,pillsC):
        print(pu,'=',pc)
