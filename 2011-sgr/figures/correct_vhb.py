import atpy, os

table = atpy.Table('combined.xml', verbose=False)

vgsr_index = table.columns.keys.index('V_{gsr} (fxcor)')


update_data = table.data
Vhb = 17.135
Vhberr = 0.165
RRLMv = +0.69

for row in update_data:
    if row[vgsr_index] > 0.:
        
        
        row[table.columns.keys.index('V_{hb}')] = Vhb
        
        row[table.columns.keys.index('V_{hb,err}')] = Vhberr
        row[table.columns.keys.index('D')] = 10**((Vhb - RRLMv  + 5 - row[table.columns.keys.index('E(B-V)')]*2.682)/5.) * 10**(-3)

        row[table.columns.keys.index('D_{err}')] = 10**((Vhb + Vhberr - RRLMv + 5 - row[table.columns.keys.index('E(B-V)')]*2.682)/5.) * 10**(-3)
        row[table.columns.keys.index('V-V_{HB}')] = row['V (Jester 2005)'] - Vhb
        row[table.columns.keys.index('[Fe/H] (Battaglia)')] = -2.81 + 0.44*(row['EW_{CaT2}'] + row['EW_{CaT3}'] + 0.64*(row['V (Jester 2005)'] - Vhb))
        row[table.columns.keys.index('[Fe/H]_{err} (Battaglia)')] = 0.16 + 0.04*(0.02*(row['V (Jester 2005)'] - (Vhb + Vhberr)))
        
        
        
        

        
        
os.system('rm -f combined.xml') 
table.write('combined.xml', verbose=False)