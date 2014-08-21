# -*- coding: cp1252 -*-
import os, sys, glob, stat, getopt, platform, subprocess, multiprocessing, fnmatch, re
import Tkinter, Tkconstants, tkFileDialog, time, codecs, xml.etree.ElementTree as ET
import shutil

cziStr          = '.czi'
cziFilterStr    = '*.czi'
xformFileStr    = '*.xform'
brightFieldStr0 = '_C0.nrrd'
Channel1Str     = '_C1.nrrd'
GdbIcpExec      = 'D:/Farsight-files/fiji-win64/Fiji.app/gdbicp.exe'
FijiExec        = 'D:/Farsight-files/fiji-win64/Fiji.app/ImageJ-win64.exe'

class TkFileDialog(Tkinter.Frame):
  def __init__(self, root):
    pform = platform.system()
    self.execExt = ''
    self.execPref = ''
    if pform.lower() == 'windows':
      self.execExt = '.exe'
    else:
      self.execPref = './'
    numCores = multiprocessing.cpu_count()
    if numCores>10:
      numCores = 10
    self.numCoresToUse = str(int( round(numCores*0.9,0) ));
    self.stdstrring = ''
    self.stdsterr   = ''

    Tkinter.Frame.__init__(self, root)

    # options for buttons
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}

    # define buttons
    Tkinter.Button(self, text='Shade Correct Set', command=self.shadecorrect).pack(**button_opt)
    Tkinter.Button(self, text='Register Set', command=self.register).pack(**button_opt)
    Tkinter.Button(self, text='Shade Correct+Register Set',
                  command=self.shadecorrectnregisterfiles).pack(**button_opt)
    Tkinter.Button(self, text='Shade Correct+Register Sets in Folders',
                  command=self.shadecorrectnregisterdirs).pack(**button_opt)

    # define options for opening or saving a file
    self.file_opt = options = {}
    options['defaultextension'] = '.czi'
    options['filetypes'] = [('czi files', '.czi'),('all files', '.*')]
    options['parent'] = root
    options['title'] = 'Select files for processing'

    # if you use the multiple file version of the module functions this option is set automatically.
    options['multiple'] = 1

  def askopenfilename( self ):
    # get filename
    root = Tkinter.Tk()
    filename = tkFileDialog.askopenfilename(**self.file_opt)
    filenamelist = root.tk.splitlist(filename)
    return filenamelist
  
  def askopendirectory( self ):
    # get filename
    root = Tkinter.Tk()
    directoryname = tkFileDialog.askdirectory()
    directorylist = root.tk.splitlist(directoryname)
    return directorylist

  def shadecorrectnregisterfiles(self):
    files = self.askopenfilename()
    print files
    self.shadecorrectnregisterloop(files)

  def shadecorrectnregisterdirs(self):
    dirs = self.askopendirectory()
    print dirs
    self.shadecorrectnregisterloop(dirs)

  def shadecorrectnregisterloop(self,files):
#    time.sleep(600)
    for filename in files:
      print filename
#      self.WriteNrrdTilesFromCziNShadeCorrect(filename)
#    class PassFileNamesFromCLInsteadOfTk():
#    def __init__(self, files):'''
    skip = 0              #No gui code
    for filename in files:
      if skip:              #
        skip = 0            #
      else:             #
        self.stdstrring = ''
        self.stdsterr   = ''
        if filename.find(cziStr)!=-1:
          self.WriteNrrdTilesFromCziNShadeCorrect(filename)
          print filename
          logfile = os.path.splitext(filename)[0]+'.log'
          errfile = os.path.splitext(filename)[0]+'.err'
          f1 = open( logfile, 'a' )
          f1.write( self.stdstrring )
          f1.close()
          f2 = open( errfile, 'a' )
          f2.write( self.stdsterr )
          f2.close()
        else:
          registerFiles = []
          for root, dirs, czifiles in os.walk(filename):
            for filenameczi in fnmatch.filter(czifiles, cziFilterStr):
              prcFilename = os.path.join(root,filenameczi)
              self.WriteNrrdTilesFromCziNShadeCorrect(prcFilename)
              logfile = os.path.splitext(prcFilename)[0]+'.log'
              errfile = os.path.splitext(prcFilename)[0]+'.err'
              f1 = open( logfile, 'a' )
              f1.write( self.stdstrring )
              f1.close()
              f2 = open( errfile, 'a' )
              f2.write( self.stdsterr )
              f2.close()
              self.stdstrring = ''
              self.stdsterr   = ''
              registerFiles.append( prcFilename )
          if len(registerFiles)>0:
            brightFeildString = 'C0_IlluminationCorrected_stitched.tif'
            channelsString = '_IlluminationCorrected_stitched.tif'
            subsampleFlag = 1
            self.registerThese( registerFiles,brightFeildString,channelsString,subsampleFlag )
            brightFeildString = 'C0_IlluminationCorrected_stitched_registered.tif'
            channelsString = '_IlluminationCorrected_stitched_registered.tif'
            subsampleFlag = 0
            self.registerThese( registerFiles,brightFeildString,channelsString,subsampleFlag )
    print('I am done processing. Close the window.')

  def registerThese( self,registerFiles,brightFeildString,channelsString,subsampleFlag ):
    #Do registration
    if len(registerFiles)<2:
      print 'Single czi passed, hence no registration required.', registerFiles
      return
    while len(registerFiles):
      #Take the first file get path
      currentPath = os.path.dirname(registerFiles[0])
      #Find all files with the same path
      currentRegList  = []
      deleteIndexList = []
      currentRegList.append( registerFiles[0] )
      for i in range( 1, len(registerFiles) ):
        if registerFiles[i].find(currentPath)!=-1:
          deleteIndexList.append( i )
          currentRegList.append( registerFiles[i] )
      #Delete the files that will be registered from the master list
      currentRegList.sort()
      deleteIndexList.append(0)
      deleteIndexList.sort()
      deleteIndexList.reverse()
      
      del registerFiles[deleteIndexList[i]]

      #Find C0 tiff from the first round of staining
      searchPattern = os.path.splitext(os.path.basename(currentRegList[0]))[0]
      brightFeildFound1 = 0
      brightFeildChannel1 = ''
      currentSubPath = os.path.dirname(currentRegList[0])
      for fileName in os.listdir(currentSubPath):
        if os.path.basename(fileName).find(searchPattern)!=-1:
          if os.path.basename(fileName).find('C0_IlluminationCorrected_stitched.tif')!=-1:
            brightFeildChannel1 = currentSubPath + '/' + fileName
            brightFeildFound1   = 1
      if brightFeildFound1:
        self.downsampleAndRunAutoContrastInImageJ( brightFeildChannel1, 0, subsampleFlag )
        if subsampleFlag:
          subsampleString = '_subsample.tif'
        else:
          subsampleString = '_cropped.tif'
        subsampleFilename1 = os.path.splitext(brightFeildChannel1)[0]+subsampleString
        for i in range( 1, len(currentRegList) ):
          #Find bright feild in subsequent rounds
          searchPattern = os.path.splitext(os.path.basename(currentRegList[i]))[0]
          brightFeildChannel2 = ''
          round2FlourChannels = []
          brightFeildFound2   = 0
          for fileName in os.listdir(currentSubPath):
            if os.path.basename(fileName).find(searchPattern)!=-1:
              if os.path.basename(fileName).find(brightFeildString)!=-1:
                brightFeildChannel2 = currentSubPath + '/' + fileName
                brightFeildFound2   = 1
              else:
                if os.path.basename(fileName).find(channelsString)!=-1:
                  fullPathToFile = currentSubPath + '/' + fileName
                  round2FlourChannels.append( fullPathToFile )
          if brightFeildFound2:
            #Start by downsampling images
            self.downsampleAndRunAutoContrastInImageJ( brightFeildChannel2, 1, subsampleFlag )
            if subsampleFlag:
              subsampleScaleFile = os.path.splitext(brightFeildChannel2)[0]+'.subsample'
              f = open(subsampleScaleFile, 'r')
              subsampleScaling = float(f.read())
              f.close()
              os.remove( subsampleScaleFile )
            else:
              subsampleScaling = 1.0
            subsampleFilename2 = os.path.splitext(brightFeildChannel2)[0]+subsampleString
            #Run gdbicp
            args = [GdbIcpExec, subsampleFilename1, subsampleFilename2, '-model', '4', '-no_render', '1']
            currentDir = os.getcwd()
            os.chdir(os.path.dirname(subsampleFilename2))
            self.RunExec( args, subsampleFilename2 )
            os.chdir(currentDir)
            os.remove( subsampleFilename2 )
            for fileName in os.listdir(currentSubPath):
              if os.path.splitext(fileName)[1].find('xform')!=-1:
                if os.path.basename(fileName).find(os.path.basename(os.path.splitext(subsampleFilename1)[0]))!=-1:
                  if os.path.basename(fileName).find(os.path.basename(os.path.splitext(subsampleFilename2)[0]))!=-1:
                    xformFile = currentSubPath + '/' + fileName
                    f = open(xformFile, 'r')
                    counter1 = 0
                    for line in f.readlines():
                      if line.find('SIMILARITY')!=-1 or counter1:
                        counter1 = counter1 + 1
                      if counter1==5:
                        stringArray = line.split()
                        xTranslate = -1*subsampleScaling*(float(stringArray[2])-float(stringArray[0]))
                        yTranslate = -1*subsampleScaling*(float(stringArray[3])-float(stringArray[1])) #Origin in ITK
                    f.close()
                    os.remove( xformFile )
                    if counter1<5:
                      print( 'Bad GDBICP file and skipping following file:' )
                      print currentRegList[i]
                    else:
                      TranslateExe = 'TranslateThenAffineRegister'+self.execExt
                      args = [ TranslateExe, brightFeildChannel1, brightFeildChannel2, str(xTranslate), str(yTranslate) ];
                      args.extend( round2FlourChannels );
                      self.RunExec( args, brightFeildChannel2 )
                    os.remove( brightFeildChannel2 )
                    for removeName in round2FlourChannels:
                      os.remove( removeName )
          else:
            print( 'Failed to find the registration brightfield and skipping following file:' )
            print currentRegList[i]
            registerFiles.clear
        os.remove( subsampleFilename1 )
      else:
        print( 'Failed to find the registration brightfield and skipping following files:' )
        print currentRegList

  #Run autocontras fiji/imagej macro
  def downsampleAndRunAutoContrastInImageJ( self, filename, writeScalingFlag, subsampleFlag ):
    if subsampleFlag:
      ResampleExe = 'ResampleImage'+self.execExt
    else:
      ResampleExe = 'CropForReg'+self.execExt
    if writeScalingFlag and subsampleFlag:
      args = [ ResampleExe, filename, '1' ]
    else:
      args = [ ResampleExe, filename ]
    self.RunExec( args, filename )
    subsampleFilename = os.path.splitext(filename)[0]+'_subsample.tif'
    ijFile = os.path.splitext(os.path.basename(filename))[0]+'.ijm'
    ijmacroFile =  os.path.dirname(FijiExec)+'/macros/'+ijFile #Fiji is defaulting to FIJI_EXE_PATH/macros/
    f = open(ijmacroFile, 'w')
    f.write('open("')
    subsampleFilenameIJ = subsampleFilename.replace('/', '\\\\')
    f.write(subsampleFilenameIJ)
    f.write('");')
    f.write('\nrun("Enhance Contrast", "saturated=0.35 normalize");\nrun("Save");\nclose();')
    f.close()
    #Run imagej/fiji to enhance contrast
    args = [FijiExec,'--headless','-macro',ijFile] #Fiji is defaulting to FIJI_EXE_PATH/macros/
    self.RunExec( args, ijmacroFile )
    os.remove( ijmacroFile )

  def register( self ):
    #files = self.askopenfilename()
    #for filename in files:
    #    print filename 
    print('Nothing implemented yet. Empty for now.')

  def shadecorrect( self ):
    #files = self.askopenfilename()
    #for filename in files:
    #    print filename 
    print('Nothing implemented yet. Empty for now.')

  def WriteNrrdTilesFromCziNShadeCorrect( self, filename ):
    #Get Metadata
    self.RunBioformatsMetaReader(filename)
    #Nrrd files before conversion
    searchStr = os.path.splitext(filename)[0] + '*.nrrd'
    #Write Nrrd files
    self.findChannelTifsAndWriteNrrds( filename )
    #Search and stitch new nrrd files into dicoms
    IlluminationEx = self.execPref+'Illumination'+self.execExt
    dicomConverter = self.execPref+'NrrdToDicom'+self.execExt
    filesThatExistPost = glob.glob(searchStr)
	rowLength = getRowLength( filename )
    for nrrdFile in filesThatExistPost:
      if nrrdFile.find(brightFieldStr0)!= -1:
        args = [ IlluminationEx, nrrdFile, '0', self.numCoresToUse, '1' ]
      else:
        args = [ IlluminationEx, nrrdFile, rowLength, self.numCoresToUse ]
##      args1 = [ dicomConverter, nrrdFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
##      if nrrdFile.find(brightFieldStr0)!=-1:
##        self.RunExec( args1, nrrdFile ) #stitches raw uncorrected images
      self.RunExec( args, nrrdFile )
      nrrdIllFile = os.path.splitext(nrrdFile)[0]+'_IlluminationCorrected.nrrd'
      if os.path.exists(nrrdIllFile):
        args = [ dicomConverter, nrrdIllFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
        self.RunExec( args, nrrdIllFile )
        os.remove( nrrdFile )
        if os.path.exists(os.path.splitext(nrrdIllFile)[0]+'_stitched.tif')==0:
          self.stdsterr += 'Stitching failed on file: '+nrrdIllFile+'\n'
        else:
          os.remove( nrrdIllFile )
      else:
        self.stdsterr   += 'Illumination correction failed on file: '+nrrdFile+'\n'

  def RunBioformatsMetaReader( self, filename ):
    xmlFilename = os.path.splitext(filename)[0]+'.xml'
    errFilename = os.path.splitext(filename)[0]+'.errlog'
    javaExec = 'C:\\Program Files\\Java\\jre8\\bin\\java.exe' #'java'+self.execExt
    args = [ javaExec , '-mx512m', '-cp', 'loci_tools.jar', 'loci.formats.tools.ImageInfo',
             '-nopix', '-omexml-only', filename ]
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    self.stdsterr   += stderr
    f = open( xmlFilename, 'w' )
    f.write(stdout)
    f.close()

  def RunExec( self, args, filename ):
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    self.stdstrring += stdout
    self.stdsterr   += stderr

  def findChannelTifsAndWriteNrrds( self, filename ):
    channelNamesArray = ['c01','c02','c03','c04','c05','c06','c07','c08','c09','c10','c11']
    outputNamesArray  = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10']
    allFilesInDir = os.listdir(os.path.dirname(filename))
    for i in range( 0, len(channelNamesArray) ):
      channelSearchStr = os.path.splitext(os.path.basename(filename))[0] + '*' + channelNamesArray[i] + '*.tif'
      currentChannelFiles = fnmatch.filter( allFilesInDir, channelSearchStr )
      currentChannelFiles.sort()
      channelListFile = os.path.splitext(filename)[0]+'_list_'+outputNamesArray[i]+'.txt'
      f = open(channelListFile, 'w')
      for j in range( 0, len(currentChannelFiles) ):
        loopFileName = os.path.dirname(filename)+'/'+currentChannelFiles[j]
        loopFileName = loopFileName.replace('/', '\\\\')
        f.write(loopFileName)
        f.write('\n')
      f.close()
      ijFile = os.path.splitext(os.path.basename(filename))[0]+channelNamesArray[i]+'.ijm'
      ijmacroFile =  os.path.dirname(FijiExec)+'/macros/'+ijFile #Fiji is defaulting to FIJI_EXE_PATH/macros/
      f = open(ijmacroFile, 'w')
      f.write('run("Stack From List...", "open=')
      channelListFile = channelListFile.replace('\\', '\\\\')
      channelListFile = channelListFile.replace('/', '\\\\')
      f.write(channelListFile)
      f.write('");\nrun("Nrrd ... ", "nrrd=')
      channelOutputName = os.path.splitext(filename)[0]+'_'+outputNamesArray[i]+'.nrrd'
      channelOutputName = channelOutputName.replace('\\', '\\\\')
      channelOutputName = channelOutputName.replace('/', '\\\\')
      f.write(channelOutputName)
      f.write('");\n close();\nrun("Quit");')
      f.close()
      #Run imagej/fiji to enhance contrast
      args = [FijiExec,'-macro',ijFile] #,'--headless'Fiji is defaulting to FIJI_EXE_PATH/macros/
      self.RunExec( args, ijmacroFile )
      os.remove( ijmacroFile )
      os.remove( channelListFile )
      for j in range( 0, len(currentChannelFiles) ):
        loopFileName = os.path.dirname(filename)+'/'+currentChannelFiles[j]
        os.remove(loopFileName)

  def getRowLength( self, filename ):
    xmlFile = os.path.splitext(filename)[0]+'.xml'
    #First fix special charecters
    f1 = open(xmlFile)
    fileData = f1.read()
    f1.close()
    fileData = re.sub('µ','mu', fileData)
    fileData = re.sub('²','^2', fileData)
    xmlFileTemp1 = os.path.splitext(xmlFile)[0] + '_temp1.xml'
    f1 = open(xmlFileTemp1, 'w')
    f1.write(fileData)
    f1.close()
    #Replace the first two lines
    xmlFileTemp2 = os.path.splitext(xmlFile)[0] + '_temp2.xml'
    f1 = open(xmlFileTemp1, 'r')
    f2 = open(xmlFileTemp2, 'w')
    f1.readline()
    f1.readline()
    f2.write('<OME>\n')
    shutil.copyfileobj(f1, f2)
    f1.close()
    f2.close()
    os.remove( xmlFileTemp1 )
    tree = ET.parse(xmlFileTemp2)
    root = tree.getroot()
    xset = 0;
    currentX = 0.0001;
    currentY = 0.0001;
    currentRowCount = 0
    for image in root.iter('Image'):
      for pixels in root.iter('Pixels'):
        for plane in root.iter('Plane'):
          if xset == 0:
            currentX = float(plane.get('PositionX'))
            currentY = float(plane.get('PositionY'))
            currentRowCount = currentRowCount + 1
            xset = 1
          else:
            tempX = float(plane.get('PositionX'))
            tempY = float(plane.get('PositionY'))
            if tempX == currentX and tempY == currentY:
              continue
            else:
              if tempY == currentY and tempX != currentX:
                 currentX = tempX
                 currentRowCount = currentRowCount + 1
              else:
                print ('Number of tiles per row:')
                print currentRowCount
                os.remove( xmlFileTemp2 )
                return currentRowCount
    return currentRowCount

if __name__=='__main__':
#  PassFileNamesFromCLInsteadOfTk(sys.argv)  #No gui code
  root = Tkinter.Tk()
  TkFileDialog(root).pack()
  root.mainloop()

