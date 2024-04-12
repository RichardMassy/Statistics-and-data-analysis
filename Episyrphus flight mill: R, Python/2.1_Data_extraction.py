"""
Created on Thu Jun 20 13:40:29 2019

@author: Richard Massy

"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy.signal import find_peaks

def deceleration(row): # equation representing no-load deceleration
    return row["Acceleration"] + (0.1009*np.exp(1.37*row["Speed"]))-0.1
# gap lim 0.2 -> 0.25
def check_range(index,df,extra,threshold = 0.2,gap_lim = 0.25):
    test = sum(df.loc[index:index+1,"Gap"])/(2-extra)
    comps = df.loc[df.index.isin([index-2,index-1]),"Gap"]
    if any([comp > 0.3 for comp in comps]):
            return False
    # New if
    if extra > 0 and all([comp > 2*df.loc[index,"Gap"] for comp in comps]):
        return True
    if all([abs(test-comp)/test < threshold for comp in comps]):
        return True
    return False

def check_miss(index,df,gap_lim = 0.25):
    # index + 3 and +4 changed to +2 and +3
    comps = df.loc[df.index.isin([index-2,index-1,index+2,index+3]),"Gap"]
    if (np.mean(comps) > gap_lim) | (df.loc[index,"Gap"] < 1.6*np.mean(comps)):
        return False
    # comp comparison < 0.16*c+0.005 changed to 0.2*c+0.005
    if all([all([abs(c-comp)<0.2*c+0.005 for c in comps]) for comp in comps]):
        return True
    else:
        return False

def df_filter_populate(df,folder,mill):
    """Identifying cases of erroneous opposite direction records"""
    dom_dir = df.Direction.value_counts().index[0]
    df["DCh"] = df.Direction.ne(df.Direction.shift())
    df["Gap"] = df["Time"].diff()
    sel = df[df.DCh.shift(-1) & df.DCh & (df.Direction != dom_dir)].index
    
    errors = {"dir_skip":[],"time_shift":[],"rec_miss":[],"rec_add":[],
              "unfixed":[]}
    
    """1st fix: at high speeds the first sensor is skipped = double record"""
    to_del = []
    for index in reversed(sel):
        if df.loc[index+1,"Gap"] < 0.005:
            errors["dir_skip"].append(df[index-5:index+5].copy())
            df.loc[index,"Direction"] = dom_dir
            to_del.append(index+1)
    df.drop(index = to_del,inplace=True)
    df = df.reset_index(drop=True)
    """1.5 fix: extra record but normal speed. Rec + rec+1 = desired Rec"""
    to_del = []
    sel = df[df.DCh.shift(-1) & df.DCh & (df.Direction != dom_dir)].index
    for index in reversed(sel[sel>5]):
        if check_range(index,df,extra=True):
            errors["rec_add"].append(df[index-5:index+5].copy())
            df.loc[index+1,"Gap"] = sum(df.loc[index:index+1,"Gap"])
            to_del.append(index)            
    df.drop(index = to_del,inplace=True)
    df = df.reset_index(drop=True)
            
    """3rd fix: missed timestamps, maybe coding wheel moves out of sensor"""
    sel = df[df.DCh.shift(-1) & df.DCh & (df.Direction != dom_dir)].index
    to_del = []
    to_insert = []
    for index in reversed(sel[sel>5]):
        if check_miss(index,df):
            errors["rec_miss"].append(df[index-5:index+5].copy())
            duration = sum(df.loc[index:index+1,"Gap"])
            n = round(duration/
                      np.mean(df.loc[[index-2,index-1,index+2,index+3],"Gap"]))
            insert = df.loc[[index-1]*n]
            insert["Time"]=insert["Time"]+[duration*i/n for i in range(1,n+1)]        
            insert["Gap"] = duration/n
            insert.index = np.arange(n)*0.01 + index
            to_insert.append(insert)
            to_del.extend([index,index+1,index+2])
    df.drop(index = to_del,inplace=True)
    df = pd.concat([df]+to_insert).sort_index()
    df = df.reset_index(drop=True)

    """Elimination of peaks (check_range without DCh), 1st pass: raw speed"""
    peaks,_ = find_peaks(0.03534292 / df["Gap"],threshold=0.2)
    peaksInfo = {"peaks1N":len(peaks)}
    df.loc[peaks,"treat"] = df.Gap.rolling(2).mean()
    df.loc[peaks+1,"treat"] = df.Gap.shift().rolling(2).mean()
    df["Gap"] = np.where(BA:=(((df.treat-df.Gap.shift(2))/df.treat)<0.2),
                         df.treat,df.Gap)
    peaksInfo["peaks1T"] = sum(BA)
    
    """2nd fix: time is shifted between recordings, mean = normal time"""
    sel = df[df.DCh.shift(-1) & df.DCh & (df.Direction != dom_dir)].index
    for index in reversed(sel[sel>5]):
        if check_range(index,df,extra=False):
            errors["time_shift"].append(df[index-5:index+5].copy())
            to_insert = {"Time":np.mean(df.loc[[index-1,index+1],"Time"]),
                         "Gap":np.mean(df.loc[index:index+1,"Gap"])}
            df.loc[index,["Time","Direction"]] = [to_insert["Time"],dom_dir]
            df.loc[index:index+1,"Gap"] = to_insert["Gap"]
        else:
            errors["unfixed"].append(df[index-5:index+5].copy())

    """2nd pass: find peaks on residual speed (remainder after smoothing)"""
    df["SpeedR"] = 0.03534292 / df["Gap"]
    df["resM"] = abs(df.SpeedR - df.SpeedR.rolling(16,center=True).mean())
    peaks,_ = find_peaks(df["resM"],prominence=0.08,width=(0,3))
    peaksInfo["peaks2N"] = len(peaks)
    df.loc[peaks,"treat2"] = df.Gap.rolling(2).mean()
    df.loc[peaks-1,"treat2"] = df.Gap.shift(-1).rolling(2).mean()
    df["Gap"] = np.where(BA:=(((df.treat2-df.Gap.shift(2))/df.treat2<0.2)&
                              (df.Gap.shift(2)<0.2)),df.treat2,df.Gap)
    peaksInfo["peaks2T"] = sum(BA)

    """Populating df"""
    df["SpeedR"] = 0.03534292 / df["Gap"]
    df["Speed"] = df["SpeedR"].rolling(2).mean()
    df["Acceleration"] = df["Speed"].diff() / df["Gap"]
    df["A_insect"] = df.apply(deceleration,axis=1)
    df["residual"] = df.Speed - df.Speed.rolling(16,center=True).mean()
    df["non_repeat"] = df.residual.diff(16)
    df["Roll_Avg_Speed"] = df.Speed.rolling(6,center=True).mean().dropna()
    """Identify flights: remove rocking and slow movement between flights"""
    df["flying"] = (df.Speed.ge(0.1) & (df.Direction == dom_dir) &
                    ((df.A_insect > 0) | df.Speed.ge(0.2)))
    df["flight"] = (df["flying"]==False).cumsum().where(df["flying"]==True)
    df = df[3:]
    
    """Aggregate statistics on separate flight events"""
    groups = df.groupby("flight")
    flights = groups.Speed.describe().add_prefix("Speed_")
    flights.rename(columns={"Speed_count":"Count"},inplace=True)
    
    flights = flights.join(groups.agg(
        Distance=("Time",lambda x: len(x)*0.03534292),
        Duration=("Time",lambda x: max(x)-min(x)),
        Accel_mean=("Acceleration","mean"),
        A_insect_mean=("A_insect","mean"),
        Residual_mov=("residual",sum),
        T_start=("Time",min),
        T_end=("Time",max)))

    def flights_2nd_stats(x):
        res_tot = abs(x.residual).sum(skipna=True)
        non_rep_tot = abs(x.non_repeat).sum(skipna=True)
        return pd.Series([np.ma.average(x["Speed"],weights=x["Gap"]),
                          np.ma.average(x["Acceleration"],weights=x["Gap"]),
                          np.ma.average(x["A_insect"],weights=x["Gap"]),
                          res_tot/non_rep_tot,(res_tot-non_rep_tot)/len(x)],
                         index=["Speed","Acceleration","A_insect","Behav_rep",
                                "Behav_add"])
    flights = flights.join(groups.apply(flights_2nd_stats))
    flights = flights.assign(Fly = folder.stem+mill,dom_dir = dom_dir)
    flights = flights[flights.Count.ge(64) & flights.Duration.ge(10)] # Filter
    
    """Stats across all movement"""
    d_counts = df.Direction.value_counts()
    stats = {"Fly":folder.stem+mill,
             "Distance":flights.Distance.sum(),
             "Acceleration":np.ma.average(flights["Acceleration"],
                                          weights=flights["Duration"]),
             "A_insect":np.ma.average(flights["A_insect"],
                                      weights=flights["Duration"]),
             "Max_speed":flights.Speed_max.max(),
             "Flight_time":flights.Duration.sum(),
             "Longest_flight":flights.Duration.max(),
             "N_flights":len(flights),
             "Last_flight_end":flights.T_end.max(),
             "F_counts":d_counts.get("F",0),
             "R_counts":d_counts.get("R",0),
             "dom_dir":dom_dir,
             "Residual_mov":flights.Residual_mov.sum,
             "Behav_rep":np.ma.average(flights["Behav_rep"],
                                       weights=flights["Duration"]),    
             "Behav_add":np.ma.average(flights["Behav_add"],
                                       weights=flights["Duration"])
             }|{"DCh_{}".format(name):len(er) for name,er in errors.items()
                }|peaksInfo
    
    stats.update({"Percent_flying":stats["Flight_time"] / 144,
                  "Mean_speed":stats["Distance"] / stats["Flight_time"],
                  "Total_dist":abs(stats["F_counts"]-stats["R_counts"]),
                  "FR_ratio":max(stats["F_counts"]/stats["R_counts"],
                                 stats["R_counts"]/stats["F_counts"])})
    
    """Experiment time transect"""    
    # Changed from <1 to >= 59
    df["Transect"] = np.where((df.Time%60 >=59) & df.flight.notna(),
                              np.rint(df.Time/60),np.nan)
    df["Transect_nf"] = np.where((df.Time%60 >=59) & df.flight.isna(),
                                 np.rint(df.Time/60),np.nan)
    transect = df.Transect.value_counts()
    transect_nf = df.Transect_nf.value_counts()
    transect = transect.drop(transect.index.intersection(transect_nf.index))
    
    """Transform errors from dictionary of lists of dfs to df"""
    Errors=pd.concat([pd.DataFrame()]+[pd.concat(dfs).assign(Error = name)
                                       for name,dfs in errors.items() if dfs])
    
    return flights, stats, Errors, transect

"""NEW 2.62"""
def mill_analysis(expDirs,home=""):
    for expDir in [expDirs] if isinstance(expDirs,str) else expDirs:
        expDir = Path(home+expDir)
        (expDir/"Old").mkdir(exist_ok=True)
        for file in ["All flights.csv","Flight_stats.csv","Errors.csv",
                     "Transects.csv"]:
            if (expDir/("Old/"+file)).exists():
                (expDir/("Old/"+file)).replace(expDir/("Old/"+file+".old"))
            if (expDir/file).exists():
                (expDir/file).replace(expDir/("Old/"+file))
        flights = []
        stats = []
        errors = {}
        transects = {}
        for folder in (expDir/"Mills").glob("[0-9]"*4):
            print(folder)
            for mill in ["mill1","mill2","mill3","mill4"]:
                if not (files:=list(folder.glob("*{}(*)*.csv".format(mill)))):
                    continue
                file = files[np.argmax([file.stem[-2] for file in files])]
                
                df = pd.read_csv(file,header=0,skipfooter=1,engine="python")
                
                if len(df) < 1600:
                    print(file.stem+" <1600 data points")
                    continue
                print(file.stem)
                
                output = df_filter_populate(df,folder=folder,mill=mill)
                
                flights.append(output[0])
                stats.append(output[1])
                errors[folder.stem+mill] = output[2]
                transects[folder.stem+mill] = output[3]
    
        Flights = pd.concat(flights)
        Stats = pd.DataFrame(stats)
        Errors = pd.concat([pd.DataFrame()]+[df.assign(Fly=name) for name, df 
                            in errors.items() if not df.empty])
        Transects=pd.DataFrame.from_dict(transects).reindex(range(240),
                                                            fill_value=np.nan)
        Flights.to_csv(expDir/"All flights.csv")
        Stats.to_csv(expDir/"Flight_stats.csv",index=False)
        Errors.to_csv(expDir/"Errors.csv",index_label="Index")
        Transects.to_csv(expDir/"Transects.csv",index_label="Index")
        
        print("Finished {}".format(expDir.stem))
        
#%%
"Run the extraction process"
import timeit
starting_time = timeit.default_timer()
mill_analysis(["21 Autumn","22 Summer"],home="~/Data S2")
print("Time difference :", timeit.default_timer() - starting_time)
"""
End
of
Script
"""