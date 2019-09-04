'''
Created on 28 Oct 2018

@author: thomasgumbricht
'''
import os
import numpy as np
from geoimagine.mask import MultiBandMasking, SingleBandMasking
#gc = grabarge collect
import gc

class ProcessOverlay:
    '''class for overlay processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        if self.process.proc.processid in ['subtractseasonsancillary']:
            self._GetLayerMask()

        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                dstCompL = []
                contFlag = False
                for dstcomp in self.process.dstLayerD[locus][datum]:
                    dstCompL.append(dstcomp)
                    if not self.process.dstLayerD[locus][datum][dstcomp]._Exists() or self.process.overwrite:
                        pass
                    else:
                        if self.verbose > 1:
                            print ('registering',self.process.dstLayerD[locus][datum][dstcomp].FPN)
                        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)
                        
                        #The locus and datum is done, clean the class variables, and set contFlag to True
                        self.process.dstLayerD[locus][datum][dstcomp] = None
                        for srccomp in self.process.srcLayerD[locus][datum]:
                            self.process.srcLayerD[locus][datum][srccomp] = None
                        #grabage collect
                        gc.collect()
                        if self.verbose > 1:
                            print ('    done, continuing',locus, datum, dstCompL[0])
                        contFlag = True
                if contFlag:
                    continue
                srcCompL = []
                missing = False 
                layeridD = {}
                for srccomp in self.process.srcLayerD[locus][datum]:
                    if not (self.process.srcLayerD[locus][datum][srccomp]):
                        #The locus and datum is missinbg, clean the class variables, and set contFlag to True
                        print ('composition missing',datum)
                        for srccomp in self.process.srcLayerD[locus][datum]:
                            self.process.srcLayerD[locus][datum][srccomp] = None
                        for dstcomp in self.process.dstLayerD[locus][datum]:
                            self.process.dstLayerD[locus][datum][dstcomp] = None
                        #grabage collect
                        gc.collect()
                        contFlag = True
                    else: 
                        if os.path.exists(self.process.srcLayerD[locus][datum][srccomp].FPN):
                            layerId = self.process.srcLayerD[locus][datum][srccomp].comp.id
                            layeridD[layerId] = srccomp
                            self.process.srcLayerD[locus][datum][srccomp].ReadRasterLayer()
                            srcCompL.append(srccomp)
                        else:
                            print ('missing file',self.process.srcLayerD[locus][datum][srccomp].FPN)    
                            BALLE
                if missing:
                    #For conditional overlay, just duplicate the existing layer 
                    if len(srcCompL) > 0:
                        for srccomp in self.process.srcLayerD[locus][datum]:
                            if not srccomp in srcCompL:
                                pass
                            SNULLEBULLE
                    else:
                        print ('no input for date',datum)
                        contFlag = True
                if contFlag:
                    continue
                print ('    processing',locus, datum, dstCompL[0])
                if self.process.proc.processid in ['conditionalsmapoverlay']:  
                    #Conditional overlay expects 1 or 2 source layers and 1 dst layer
                    self._ConditionalOverlay(locus, datum, srcCompL, dstCompL[0]) 
                elif 'fractionadjust' in self.process.proc.processid:  
                    #Conditional overlay expects 1 or 2 source layers and 1 dst layer
                    self._FractionAdjust(locus, datum, srcCompL, dstCompL[0],layeridD) 
                    
                elif self.process.proc.processid in ['subtractseasonsancillary']:  
                    #Conditional overlay expects 1 or 2 source layers and 1 dst layer
                    self._SubtractSeasonsAncillary(locus, datum, dstCompL[0]) 
                elif self.process.proc.processid in ['signiftrendsancillary','signiftrendsmodis']:  
                    #Conditional overlay expects 1 or 2 source layers and 1 dst layer
                    self._SignifTrend(locus, datum, dstCompL) 
                elif self.process.proc.processid.lower()[0:14] == 'dualtrendscomp':

                    self._DualPTrendsComp(locus, datum, dstCompL[0], layeridD)
                else:
                    exitStr = 'EXITING, unknown process in ProcessOverlay (overaly.py: %(s)s)' %{'s':self.process.proc.processid}
                    exit(exitStr)
                for dstComp in dstCompL:
                    #Register the layer
                    print ('registering',self.process.dstLayerD[locus][datum][dstComp].FPN)

                    self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)
                #Clear the class instances to regain memory
                for srccomp in self.process.srcLayerD[locus][datum]:
                    self.process.srcLayerD[locus][datum][srccomp] = None
                for dstcomp in self.process.dstLayerD[locus][datum]:
                    self.process.dstLayerD[locus][datum][dstcomp] = None
                #garbage collect
                gc.collect()
                              
    def _ConditionalOverlay(self,locus, datum, srcCompL, dstComp):
        dstNull = self.process.dstLayerD[locus][datum][dstComp].comp.cellnull
        #Create single dst array
        if self.process.params.method == 'avg': 
            dstBAND = np.zeros_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #Copy array
                srcBAND = BAND*1
                #Set nan to zero
                srcBAND[np.isnan(srcBAND)] = 0
                dstBAND += srcBAND   
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull
            #Set to average where both input bands are valid
            dstBAND[ (~np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (~np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] /= 2
        
        elif self.process.params.method == 'min':
            dstBAND = np.full_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND, np.nan )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #fmin
                dstBAND = np.fmin(dstBAND, BAND )
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull
             
        elif self.process.params.method == 'max':
            dstBAND = np.full_like( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND, np.nan )
            #Loop over the scr compositions
            for c in srcCompL:
                #Set null to nan
                BAND = self.process.srcLayerD[locus][datum][c].layer.NPBAND
                cellnull = self.process.srcLayerD[locus][datum][c].comp.cellnull
                BAND[BAND == cellnull] = np.nan
                #fmax
                dstBAND = np.fmax(dstBAND, BAND )
            #Set to nan where both input layers are non
            dstBAND[ (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[0]].layer.NPBAND ) ) & (np.isnan( self.process.srcLayerD[locus][datum][srcCompL[1]].layer.NPBAND ) ) ] = dstNull

        else:   
            ERRORIGEN

        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstBAND
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][c].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()
        #Register the layer
        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstComp], self.process.overwrite, self.process.delete)
        
    def _FractionAdjust(self,locus, datum, srcCompL, dstComp,layeridD):
        
        srcBAND = self.process.srcLayerD[locus][datum][layeridD['layer1']].layer.NPBAND
        FRAC = self.process.srcLayerD[locus][datum][layeridD['fraction']].layer.NPBAND

        if self.process.params.nullaszero:
            FRAC[FRAC == self.process.srcLayerD[locus][datum][layeridD['fraction']].layer.cellnull] = 0

        if self.process.params.method[0:10] == 'insidefrac':
            dstBAND = srcBAND * FRAC
        elif self.process.params.method[0:11] == 'outsidefrac':
            dstBAND = srcBAND * (1-FRAC)
        elif self.process.params.method[0:17] == 'collectinsidefrac':
            dstBAND = srcBAND * (1/FRAC)
            mask = (FRAC == 1)
            dstBAND[mask] = srcBAND[mask]
        elif self.process.params.method[0:18] == 'collectoutsidefrac':
            dstBAND = srcBAND * (1/(1-FRAC))
            mask = (FRAC == 1)
            dstBAND[mask] = srcBAND[mask]
        if self.process.params.dstmax > -32768:
            dstBAND[dstBAND > self.process.params.dstmax] = self.process.params.dstmax
        if self.process.params.dstweight < 1 and self.process.params.dstweight > 0:
            dstBAND = dstBAND*self.process.params.dstweight + srcBAND*(1-self.process.params.dstweight)
            
            
        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstBAND
        #Mask the dst BAND
        if self.process.params.nullaszero:
            SingleBandMasking(self.process.srcLayerD[locus][datum][layeridD['layer1']], self.process.dstLayerD[locus][datum][dstComp])
        else:
            MultiBandMasking(self.process.srcLayerD[locus][datum], self.process.dstLayerD[locus][datum])
        
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][layeridD['layer1']].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()


        
    def _SubtractSeasonsAncillary(self, locus, datum, dstComp):
        '''
        '''
        srcBAND = self.process.srcLayerD[locus][datum][self.layercomp].layer.NPBAND
        SEASON = self.process.srcLayerD[locus][datum][self.seasoncomp].layer.NPBAND
        
        if self.process.params.balance == 'diff': 
            
            dstBAND = srcBAND - SEASON

        elif self.process.params.balance == 'surplus':
            dstBAND = srcBAND - SEASON
            
            #remove any negative numbers
            dstBAND *= (dstBAND>0) #Should be the fastest
                         
        elif self.process.params.balance == 'deficit':
            dstBAND = (np.minimum(srcBAND, SEASON) - SEASON)

        else:   
            ERRORIGEN
                
        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = dstBAND
        #Mask the dst BAND

        MultiBandMasking(self.process.srcLayerD[locus][datum], self.process.dstLayerD[locus][datum])
        
        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][self.layercomp].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()
        
    def _SignifTrend(self,locus,datum,dstCompL):
        '''
        '''
        if self.process.srcperiod.timestep in ['timespan-A','timespan-AS']:
            self.periods = self.process.srcperiod.enddate.year-self.process.srcperiod.startdate.year

        BAND = {}

        BANDsignif = self.process.srcLayerD[locus][datum][self.process.srcIdDict['significance']].layer.NPBAND
        BANDslope = self.process.srcLayerD[locus][datum][self.process.srcIdDict['slope']].layer.NPBAND
        BANDintercept = self.process.srcLayerD[locus][datum][self.process.srcIdDict['intercept']].layer.NPBAND
        
        BAND['insignificant'] = np.ones_like(BANDsignif)
                
        threshold = self.process.params.threshold
        if self.process.params.neg and self.process.params.pos:
            BAND['insignificant'][ (BANDsignif > -threshold) & (BANDsignif < threshold) ] = 0
        elif self.process.params.neg:
            BAND['insignificant'][ (BANDsignif > -threshold)] = 0
        elif self.process.params.pos:
            BAND['insignificant'][ (BANDsignif < threshold) ] = 0
            
        BAND['slopep'] = BANDslope*BAND['insignificant']
        
        BAND['initial'] = BANDintercept*1
        
        BAND['initialp'] = BANDintercept*BAND['insignificant']
        

        BAND['final'] = BAND['initial']+BANDslope*self.periods #Calculate change from first to last
        
        BAND['finalp'] = (BAND['initial']+BANDslope*self.periods)*BAND['insignificant'] #Calculate change from first to last
        
        BAND['change'] = BAND['final']-BAND['initial']
        BAND['changep'] = BAND['finalp']-BAND['initialp']

        for dstComp in dstCompL:
            
            #Create the dst layer
            self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
            #Set the np array as the band
            self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = BAND[dstComp]
            
        #Mask the dst BAND
        MultiBandMasking(self.process.srcLayerD[locus][datum], self.process.dstLayerD[locus][datum])
            
        for dstComp in dstCompL:
            #copy the geoformat from the src layer
            self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][self.process.srcIdDict['slope']].layer)
            #write the results
            self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray()    
    
    def _DualPTrendsComp(self, locus, datum, dstComp, layeridD):
        '''
        '''
        BANDp1 = self.process.srcLayerD[locus][datum][layeridD['ptrend1']].layer.NPBAND
        BANDp2 = self.process.srcLayerD[locus][datum][layeridD['ptrend2']].layer.NPBAND
        p1Cellnull = self.process.srcLayerD[locus][datum][layeridD['ptrend1']].comp.cellnull
        p2Cellnull = self.process.srcLayerD[locus][datum][layeridD['ptrend2']].comp.cellnull
        dstCellnull = self.process.dstLayerD[locus][datum][dstComp].comp.cellnull
        
        BAND = np.ones(BANDp1.shape, np.uint8)
        #Set and default values to no change (=5)
        BAND *= 5
        '''
        #both negative = 1
        BAND[(BANDp1 < 0) & (BANDp2 < 0) & (BANDp1 != p1Cellnull) & (BANDp2 != p2Cellnull)] = 1
        #both positive = 9
        BAND[(BANDp1 > 0) & (BANDp2 > 0) & (BANDp1 != p1Cellnull) & (BANDp2 != p2Cellnull)] = 9
        #first trend positive second trend negative = 3
        BAND[(BANDp1 > 0) & (BANDp2 < 0) & (BANDp1 != p1Cellnull) & (BANDp2 != p2Cellnull)] = 3
        #first trend negative second trend positive = 7
        BAND[(BANDp1 < 0) & (BANDp2 > 0) & (BANDp1 != p1Cellnull) & (BANDp2 != p2Cellnull)] = 7
        '''
        #both negative = 1
        BAND[(BANDp1 < -0.1) & (BANDp2 < -0.1)] = 1
        #both positive = 9
        BAND[(BANDp1 > 0.1) & (BANDp2 > 0.1)] = 9
        #first trend positive second trend negative = 3
        BAND[(BANDp1 > 0.1) & (BANDp2 < -0.1)] = 3
        #first trend negative second trend positive = 7
        BAND[(BANDp1 < -0.1) & (BANDp2 > 0.1)] = 7
        #first trend negative second with no trend = 4
        BAND[(BANDp1 < -0.1) & (BAND == 5)] = 4
        #first trend postive second with no trend = 6
        BAND[(BANDp1 > 0.1) & (BAND == 5)] = 6
        #second trend negative first with no trend = 2
        BAND[(BANDp2 < -0.1) & (BAND == 5)] = 2
        #second trend positive first with no trend = 8
        BAND[(BANDp2 > 0.1) & (BAND == 5)] = 8
        #nodata in either
        BAND[(BANDp1 == p1Cellnull) | (BANDp2 == p2Cellnull)] = dstCellnull
 
        #Create the dst layer
        self.process.dstLayerD[locus][datum][dstComp].layer = lambda:None
        #Set the np array as the band
        self.process.dstLayerD[locus][datum][dstComp].layer.NPBAND = BAND

        #copy the geoformat from the src layer
        self.process.dstLayerD[locus][datum][dstComp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][layeridD['ptrend1']].layer)
        #write the results
        self.process.dstLayerD[locus][datum][dstComp].CreateDSWriteRasterArray() 

    def _GetLayerMask(self): 
        '''Identifies the bands representing layer and mask
        '''
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                for srccomp in self.process.srcLayerD[locus][datum]:
                    print (locus,datum,srccomp)
                    if self.process.srcLayerD[locus][datum][srccomp].comp.id == 'layer':
                        self.layercomp = srccomp
                    elif self.process.srcLayerD[locus][datum][srccomp].comp.id == 'season':
                        self.seasoncomp = srccomp
                    else:
                        print (self.process.srcLayerD[locus][datum][srccomp].comp.id)
                        ERRORIGEN
                return