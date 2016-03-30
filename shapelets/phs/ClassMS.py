import numpy as np
from casacore.tables import table
from rad2hmsdms import rad2hmsdms
import ModColor
import reformat

class ClassMS():
    def __init__(self,MSname,Col="DATA",zero_flag=True,ReOrder=False,EqualizeFlag=False,DoPrint=True,DoReadData=True,
                 TimeChunkSize=None,GetBeam=False,RejectAutoCorr=False,SelectSPW=None,DelStationList=None):


        if MSname=="": exit()
        MSname=reformat.reformat(MSname,LastSlash=False)
        self.MSName=MSname
        self.ColName=Col
        self.zero_flag=zero_flag
        self.ReOrder=ReOrder
        self.EqualizeFlag=EqualizeFlag
        self.DoPrint=DoPrint
        self.TimeChunkSize=TimeChunkSize
        self.RejectAutoCorr=RejectAutoCorr
        self.SelectSPW=SelectSPW
        self.DelStationList=DelStationList
        self.ReadMSInfo(MSname,DoPrint=DoPrint)

        if DoReadData: self.ReadData()
        #self.RemoveStation()

        self.SR=None
        if GetBeam:
            self.LoadSR()




    def LoadSR(self):
        if self.SR!=None: return
        import lofar.stationresponse as lsr
        f=self.ChanFreq.flatten()
        if f.shape[0]>1:
            t=table(self.MSName+"/SPECTRAL_WINDOW/",readonly=False)
            c=t.getcol("CHAN_WIDTH")
            c.fill(np.abs((f[0:-1]-f[1::])[0]))
            t.putcol("CHAN_WIDTH",c)
            t.close()
        self.SR = lsr.stationresponse(self.MSName)
        self.SR.setDirection(self.rarad,self.decrad)
        

    def GiveBeam(self,time,ra,dec):
        Beam=np.zeros((ra.shape[0],self.na,1,2,2),dtype=np.complex)
        for i in range(ra.shape[0]):
            self.SR.setDirection(ra[i],dec[i])
            Beam[i]=self.SR.evaluate(time)
        return Beam

    def GiveMappingAnt(self,ListStrSel,(row0,row1)=(None,None),FlagAutoCorr=True):
        #ListStrSel=["RT9-RTA", "RTA-RTB", "RTC-RTD", "RT6-RT7", "RT5-RT*"]

        print ModColor.Str("  ... Building BL-mapping for %s"%str(ListStrSel))

        if row1==None:
            row0=0
            row1=self.nbl
        A0=self.A0[row0:row1]
        A1=self.A1[row0:row1]
        MapOut=np.ones((self.nbl,),dtype=np.bool)
        if FlagAutoCorr:
            ind=np.where(A0==A1)[0]
            MapOut[ind]=False


        for blsel in ListStrSel:
            if "-" in blsel:
                StrA0,StrA1=blsel.split("-")
                NumA0=np.where(np.array(self.StationNames)==StrA0)[0]
                NumA1=np.where(np.array(self.StationNames)==StrA1)[0]
                C0=((A0==NumA0)&(A1==NumA1))
                C1=((A1==NumA0)&(A0==NumA1))
            else:
                NumA0=np.where(np.array(self.StationNames)==blsel)[0]
                C0=(A0==NumA0)
                C1=(A1==NumA0)
            ind=np.where(C1|C0)[0]
            MapOut[ind]=False
        self.MapSelBLs=MapOut
        return self.MapSelBLs
                



    def SelChannel(self,(start,end,step)=(None,None,None),Revert=False):
        if start!=None:
            if Revert==False:
                ind=np.arange(self.Nchan)[start:end:step]
            else:
                ind=np.array(sorted(list(set(np.arange(self.Nchan).tolist())-set(np.arange(self.Nchan)[start:end:step].tolist()))))
            self.data=self.data[:,ind,:]
            self.flag_all=self.flag_all[:,ind,:]
            shape=self.ChanFreq.shape
            self.ChanFreq=self.ChanFreq[ind]
                
        
    def ReadData(self,t0=0,t1=-1,DoPrint=False,ReadWeight=False):
        if DoPrint==True:
            print "   ... Reading MS"
        row0=0
        row1=self.F_nrows
        if t1>t0:
            t0=t0*3600.+self.F_tstart
            t1=t1*3600.+self.F_tstart
            ind0=np.argmin(np.abs(t0-self.F_times))
            ind1=np.argmin(np.abs(t1-self.F_times))
            row0=ind0*self.nbl
            row1=ind1*self.nbl

        self.ROW0=row0
        self.ROW1=row1
        self.nRowRead=row1-row0
        nRowRead=self.nRowRead

        table_all=table(self.MSName,ack=False,readonly=False)
        self.ColNames=table_all.colnames()
        SPW=table_all.getcol('DATA_DESC_ID',row0,nRowRead)
        A0=table_all.getcol('ANTENNA1',row0,nRowRead)[SPW==self.ListSPW[0]]
        A1=table_all.getcol('ANTENNA2',row0,nRowRead)[SPW==self.ListSPW[0]]
        print self.ListSPW[0]
        time_all=table_all.getcol("TIME")[SPW==self.ListSPW[0]]
        print np.max(time_all)-np.min(time_all)
        time_slots_all=np.array(sorted(list(set(time_all))))
        ntimes=time_all.shape[0]/self.nbl

        flag_all=table_all.getcol("FLAG",row0,nRowRead)[SPW==self.ListSPW[0]]
        if ReadWeight==True:
            self.Weights=table_all.getcol("WEIGHT")

        if self.EqualizeFlag:
            for i in range(self.Nchan):
                fcol=flag_all[:,i,0]|flag_all[:,i,1]|flag_all[:,i,2]|flag_all[:,i,3]
                for pol in range(4):
                    flag_all[:,i,pol]=fcol


        self.multidata=(type(self.ColName)==list)
        self.ReverseAntOrder=(np.where((A0==0)&(A1==1))[0]).shape[0]>0
        self.swapped=False

        uvw=table_all.getcol('UVW',row0,nRowRead)[SPW==self.ListSPW[0]]

        if self.ReOrder:
            vis_all=table_all.getcol(self.ColName,row0,nRowRead)
            if self.zero_flag: vis_all[flag_all==1]=0.
            if self.zero_flag: 
                noise=(np.random.randn(vis_all.shape[0],vis_all.shape[1],vis_all.shape[2])\
                           +1j*np.random.randn(vis_all.shape[0],vis_all.shape[1],vis_all.shape[2]))*1e-6
                vis_all[flag_all==1]=noise[flag_all==1]
            vis_all[np.isnan(vis_all)]=0.
            listDataSPW=[np.swapaxes(vis_all[SPW==i,:,:],0,1) for i in self.ListSPW]
            self.data=np.concatenate(listDataSPW)#np.swapaxes(np.concatenate(listDataSPW),0,1)
            listFlagSPW=[np.swapaxes(flag_all[SPW==i,:,:],0,1) for i in self.ListSPW]
            flag_all=np.concatenate(listFlagSPW)#np.swapaxes(np.concatenate(listDataSPW),0,1)
            self.uvw=uvw
            self.swapped=True
            

        else:
            self.uvw=uvw
            if self.multidata:
                self.data=[]
                for colin in self.ColName:
                    print "... read %s"%colin
                    vis_all=table_all.getcol(colin,row0,nRowRead)[SPW==self.ListSPW[0]]
                    print " shape: %s"%str(vis_all.shape)
                    if self.zero_flag: vis_all[flag_all==1]=0.
                    vis_all[np.isnan(vis_all)]=0.
                    self.data.append(vis_all)
            else:
                vis_all=table_all.getcol(self.ColName,row0,nRowRead)
                if self.zero_flag: vis_all[flag_all==1]=0.
                vis_all[np.isnan(vis_all)]=0.
                self.data=vis_all


        self.flag_all=flag_all

        if self.RejectAutoCorr:
            indGetCorrelation=np.where(A0!=A1)[0]
            A0=A0[indGetCorrelation]
            A1=A1[indGetCorrelation]
            self.uvw=self.uvw[indGetCorrelation,:]
            time_all=time_all[indGetCorrelation]
            if self.swapped:
                self.data=self.data[:,indGetCorrelation,:]
                self.flag_all=self.flag_all[:,indGetCorrelation,:]
            else:
                self.data=self.data[indGetCorrelation,:,:]
                self.flag_all=self.flag_all[indGetCorrelation,:,:]
            self.nbl=(self.na*(self.na-1))/2

        table_all.close()
        if self.DoRevertChans:
            self.data=self.data[:,::-1,:]
            self.flag_all=self.flag_all[:,::-1,:]


        self.times_all=time_all
        self.times=time_slots_all
        self.ntimes=time_slots_all.shape[0]
        self.nrows=time_all.shape[0]

        self.IndFlag=np.where(flag_all==True)
    
        #self.NPol=vis_all.shape[2]
        self.A0=A0
        self.A1=A1
        

    def SaveAllDataStruct(self):
        t=table(self.MSName,ack=False,readonly=False)

        t.putcol('ANTENNA1',self.A0)
        t.putcol('ANTENNA2',self.A1)
        t.putcol("TIME",self.times_all)
        t.putcol("TIME_CENTROID",self.times_all)
        t.putcol("UVW",self.uvw)
        t.putcol("FLAG",self.flag_all)
        for icol in range(len(self.ColName)):
            t.putcol(self.ColName[icol],self.data[icol])
        t.close()

    def RemoveStation(self):
        
        DelStationList=self.DelStationList
        if DelStationList==None: return

        StationNames=self.StationNames
        self.MapStationsKeep=np.arange(len(StationNames))
        DelNumStationList=[]
        for Station in DelStationList:
            ind=np.where(Station==np.array(StationNames))[0]
            self.MapStationsKeep[ind]=-1
            DelNumStationList.append(ind)
            indRemove=np.where((self.A0!=ind)&(self.A1!=ind))[0]
            self.A0=self.A0[indRemove]
            self.A1=self.A1[indRemove]
            self.data=self.data[indRemove,:,:]
            self.flag_all=self.flag_all[indRemove,:,:]
            self.times_all=self.times_all[indRemove,:,:]
        self.MapStationsKeep=self.MapStationsKeep[self.MapStationsKeep!=-1]
        StationNames=(np.array(StationNames)[self.MapStationsKeep]).tolist()

        na=self.MapStationsKeep.shape[0]
        self.na=na
        self.StationPos=self.StationPos[self.MapStationsKeep,:]
        self.nbl=(na*(na-1))/2+na
        

    def ReadMSInfo(self,MSname,DoPrint=True):

        ta=table(MSname+'/ANTENNA',ack=False)
        StationNames=ta.getcol('NAME')
        
        na=ta.getcol('POSITION').shape[0]
        self.StationPos=ta.getcol('POSITION')
        nbl=(na*(na-1))/2+na
        #nbl=(na*(na-1))/2
        ta.close()

        table_all=table(MSname,ack=False,readonly=False)
        SPW=table_all.getcol('DATA_DESC_ID')
        if self.SelectSPW!=None:
            self.ListSPW=self.SelectSPW
            print "dosel"
        else:
            self.ListSPW=sorted(list(set(SPW)))
        self.F_nrows=table_all.getcol("TIME").shape[0]
        F_time_all=table_all.getcol("TIME")[SPW==self.ListSPW[0]]
        #nbl=(np.where(F_time_all==F_time_all[0])[0]).shape[0]
        F_time_slots_all=np.array(sorted(list(set(F_time_all.tolist()))))
        F_ntimes=F_time_slots_all.shape[0]
        table_all.close()

        ta_spectral=table(MSname+'/SPECTRAL_WINDOW/',ack=False)
        reffreq=ta_spectral.getcol('REF_FREQUENCY')
        chan_freq=ta_spectral.getcol('CHAN_FREQ')
        self.dFreq=ta_spectral.getcol("CHAN_WIDTH").flatten()[0]
        if chan_freq.shape[0]>len(self.ListSPW):
            print ModColor.Str("  ====================== >> More SPW in headers, modifying that error....")
            chan_freq=chan_freq[np.array(self.ListSPW),:]
            reffreq=reffreq[np.array(self.ListSPW)]


        wavelength=299792458./reffreq
        NSPW=chan_freq.shape[0]
        self.ChanFreq=chan_freq
        self.Freq_Mean=np.mean(chan_freq)
        wavelength_chan=299792458./chan_freq

        if NSPW>1:
            print "Don't deal with multiple SPW yet"


        Nchan=wavelength_chan.shape[1]
        NSPWChan=NSPW*Nchan
        ta=table(MSname+'/FIELD/',ack=False)
        rarad,decrad=ta.getcol('PHASE_DIR')[0][0]
        if rarad<0.: rarad+=2.*np.pi

        radeg=rarad*180./np.pi
        decdeg=decrad*180./np.pi
        ta.close()
         
        self.DoRevertChans=False
        if Nchan>1:
            self.DoRevertChans=(self.ChanFreq.flatten()[0]>self.ChanFreq.flatten()[-1])
        if self.DoRevertChans:
            print ModColor.Str("  ====================== >> Revert Channel order!")
            wavelength_chan=wavelength_chan[0,::-1]
            self.ChanFreq=self.ChanFreq[0,::-1]

        self.na=na
        self.Nchan=Nchan
        self.NSPW=NSPW
        self.NSPWChan=NSPWChan
        self.F_tstart=F_time_all[0]
        self.F_times_all=F_time_all
        self.F_times=F_time_slots_all
        self.F_ntimes=F_time_slots_all.shape[0]
        self.dt=F_time_slots_all[1]-F_time_slots_all[0]
        self.DTs=F_time_slots_all[-1]-F_time_slots_all[0]
        self.DTh=self.DTs/3600.
        self.radec=(rarad,decrad)
        self.rarad=rarad
        self.decrad=decrad
        self.reffreq=reffreq
        self.StationNames=StationNames
        self.wavelength_chan=wavelength_chan
        self.rac=rarad
        self.decc=decrad
        self.nbl=nbl
        self.StrRA  = rad2hmsdms(self.rarad,Type="ra").replace(" ",":")
        self.StrDEC = rad2hmsdms(self.decrad,Type="dec").replace(" ",".")

        if DoPrint==True:
            print ModColor.Str(" MS PROPERTIES: ")
            print "   - File Name: %s"%ModColor.Str(self.MSName,col="green")
            print "   - Column Name: %s"%ModColor.Str(str(self.ColName),col="green")
            print "   - Pointing center: (ra, dec)=(%s, %s) "%(rad2hmsdms(self.rarad,Type="ra").replace(" ",":")\
                                                                   ,rad2hmsdms(self.decrad,Type="dec").replace(" ","."))
            print "   - Frequency = %s MHz"%str(reffreq/1e6)
            print "   - Wavelength = ",wavelength," meters"
            print "   - Time bin = %4.1f seconds"%(self.dt)
            print "   - Total Integration time = %6.2f hours"%((F_time_all[-1]-F_time_all[0])/3600.)
            print "   - Number of antenna  = ",na
            print "   - Number of baseline = ",nbl
            print "   - Number of SPW = ",NSPW
            print "   - Number of channels = ",Nchan
            print 


    def radec2lm_scalar(self,ra,dec):
        l = np.cos(dec) * np.sin(ra - self.rarad)
        m = np.sin(dec) * np.cos(self.decrad) - np.cos(dec) * np.sin(self.decrad) * np.cos(ra - self.rarad)
        return l,m

    def SaveVis(self,vis=None,Col="CORRECTED_DATA",spw=0,DoPrint=True):
        if vis==None:
            vis=self.data
        if DoPrint: print "  Writting data in column %s"%ModColor.Str(Col,col="green")
        table_all=table(self.MSName,ack=False,readonly=False)
        if self.swapped:
            visout=np.swapaxes(vis[spw*self.Nchan:(spw+1)*self.Nchan],0,1)
        else:
            visout=vis
        table_all.putcol(Col,visout,self.ROW0,self.nRowRead)
        table_all.close()
        
    def GiveUvwBL(self,a0,a1):
        vecout=self.uvw[(self.A0==a0)&(self.A1==a1),:]
        return vecout

    def GiveVisBL(self,a0,a1,col=0,pol=None):
        if self.multidata:
            vecout=self.data[col][(self.A0==a0)&(self.A1==a1),:,:]
        else:
            vecout=self.data[(self.A0==a0)&(self.A1==a1),:,:]
        if pol!=None:
            vecout=vecout[:,:,pol]
        return vecout

    def GiveVisBLChan(self,a0,a1,chan,pol=None):
        if pol==None:
            vecout=(self.data[(self.A0==a0)&(self.A1==a1),chan,0]+self.data[(self.A0==a0)&(self.A1==a1),chan,3])/2.
        else:
            vecout=self.data[(self.A0==a0)&(self.A1==a1),chan,pol]
        return vecout

    def plotBL(self,a0,a1,pol=0):
        
        import pylab
        if self.multidata:
            vis0=self.GiveVisBL(a0,a1,col=0,pol=pol)
            vis1=self.GiveVisBL(a0,a1,col=1,pol=pol)
            pylab.clf()
            pylab.subplot(2,1,1)
            #pylab.plot(vis0.real)
            pylab.plot(np.abs(vis0))
            #pylab.subplot(2,1,2)
            #pylab.plot(vis1.real)
            pylab.plot(np.abs(vis1),ls=":")
            pylab.title("%i-%i"%(a0,a1))
            #pylab.plot(vis1.real-vis0.real)
            pylab.subplot(2,1,2)
            pylab.plot(np.angle(vis0))
            pylab.plot(np.angle(vis1),ls=":")
            pylab.draw()
            pylab.show()
        else:
            pylab.clf()
            vis=self.GiveVisBL(a0,a1,col=0,pol=pol)
            pylab.subplot(2,1,1)
            #pylab.plot(np.abs(vis))
            pylab.plot(np.real(vis))
            pylab.subplot(2,1,2)
            #pylab.plot(np.angle(vis))
            pylab.plot(np.imag(vis))
            pylab.draw()
            pylab.show()

    def GiveCol(self,ColName):
        t=table(self.MSName,readonly=False,ack=False)
        col=t.getcol(ColName)
        t.close()
        return col

    def PutColInData(self,SpwChan,pol,data):
        if self.swapped:
            self.data[SpwChan,:,pol]=data
        else:
            self.data[:,SpwChan,pol]=data

    def Restore(self):
        backname="CORRECTED_DATA_BACKUP"
        backnameFlag="FLAG_BACKUP"
        t=table(self.MSName,readonly=False,ack=False)
        if backname in t.colnames():
            print "  Copying ",backname," to CORRECTED_DATA"
            #t.putcol("CORRECTED_DATA",t.getcol(backname))
            self.CopyCol(backname,"CORRECTED_DATA")
            print "  Copying ",backnameFlag," to FLAG"
            self.CopyCol(backnameFlag,"FLAG")
            #t.putcol(,t.getcol(backnameFlag))
        t.close()

    def ZeroFlagSave(self,spw=0):
        self.flag_all.fill(0)
        if self.swapped:
            flagout=np.swapaxes(self.flag_all[spw*self.Nchan:(spw+1)*self.Nchan],0,1)
        else:
            flagout=self.flag_all
        t=table(self.MSName,readonly=False,ack=False)
        t.putcol("FLAG",flagout)
        
        t.close()

    def CopyCol(self,Colin,Colout):
        t=table(self.MSName,readonly=False,ack=False)
        if self.TimeChunkSize==None:
            print "  ... Copying column %s to %s"%(Colin,Colout)
            t.putcol(Colout,t.getcol(Colin))
        else:
            print "  ... Copying column %s to %s"%(Colin,Colout)
            TimesInt=np.arange(0,self.DTh,self.TimeChunkSize).tolist()
            if not(self.DTh in TimesInt): TimesInt.append(self.DTh)
            for i in range(len(TimesInt)-1):
                t0,t1=TimesInt[i],TimesInt[i+1]
                print t0,t1
                t0=t0*3600.+self.F_tstart
                t1=t1*3600.+self.F_tstart
                ind0=np.argmin(np.abs(t0-self.F_times))
                ind1=np.argmin(np.abs(t1-self.F_times))
                row0=ind0*self.nbl
                row1=ind1*self.nbl
                NRow=row1-row0
                t.putcol(Colout,t.getcol(Colin,row0,NRow),row0,NRow)
        t.close()

        
    def PutBackupCol(self,back="CORRECTED_DATA"):
        backname="%s_BACKUP"%back
        backnameFlag="FLAG_BACKUP"
        self.PutCasaCols()
        t=table(self.MSName,readonly=False,ack=False)
        JustAdded=False
        if not(backname in t.colnames()):
            print "  Putting column ",backname," in MS"
            desc=t.getcoldesc("CORRECTED_DATA")
            desc["name"]=backname
            desc['comment']=desc['comment'].replace(" ","_")
            t.addcols(desc)
            print "  Copying CORRECTED_DATA in CORRECTED_DATA_BACKUP"
            self.CopyCol("CORRECTED_DATA",backname)
            #t.putcol(backname,t.getcol("CORRECTED_DATA"))
        if not(backnameFlag in t.colnames()):
            desc=t.getcoldesc("FLAG")
            desc["name"]=backnameFlag
            desc['comment']=desc['comment'].replace(" ","_")
            t.addcols(desc)
            self.CopyCol("FLAG",backnameFlag)
            #t.putcol(backnameFlag,t.getcol("FLAG"))
            JustAdded=True
        #else:
            #print "  Column %s already there..."%backname
        t.close()
        return JustAdded

    def PutNewCol(self,Name,LikeCol="CORRECTED_DATA"):
        if not(Name in self.ColNames):
            print "  Putting column %s in MS, with format of %s"%(Name,LikeCol)
            t=table(self.MSName,readonly=False,ack=False)
            desc=t.getcoldesc(LikeCol)
            desc["name"]=Name
            t.addcols(desc) 
            t.close()
    
    def RotateMS(self,radec):
        import ModRotate
        ModRotate.Rotate(self,radec)
        ta=table(self.MSName+'/FIELD/',ack=False,readonly=False)
        ra,dec=radec
        radec=np.array([[[ra,dec]]])
        ta.putcol("DELAY_DIR",radec)
        ta.putcol("PHASE_DIR",radec)
        ta.putcol("REFERENCE_DIR",radec)
        ta.close()
        t=table(self.MSName,ack=False,readonly=False)
        t.putcol(self.ColName,self.data)
        t.putcol("UVW",self.uvw)
        t.close()
    
    def PutCasaCols(self):
        import casacore.tables
        casacore.tables.addImagingColumns(self.MSName,ack=False)
        #self.PutNewCol("CORRECTED_DATA")
        #self.PutNewCol("MODEL_DATA")
