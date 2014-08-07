import os, sys, glob, stat, getopt, platform, subprocess, multiprocessing, fnmatch
import Tkinter, Tkconstants, tkFileDialog, time

cziStr          = '.czi'
cziFilterStr    = '*.czi'
xformFileStr    = '*.xform'
brightFieldStr0 = '_C0.nrrd'
Channel1Str     = '_C1.nrrd'
Channel3Str     = '_C3.nrrd'
GdbIcpExec      = 'D:/fiji/gdbicp.exe'
FijiExec        = 'D:/fiji/ImageJ-win64.exe'

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
          #Do registration
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
            for i in range( 0, len(deleteIndexList) ):
              del registerFiles[deleteIndexList[i]]

            #Find C0 tiff from the first round of staining
            searchPattern = os.path.splitext(os.path.basename(currentRegList[0]))[0]
            brightFeildString = 'C0_IlluminationCorrected_stitched.tif'
            brightFeildFound1 = 0
            brightFeildChannel1 = ''
            currentSubPath = os.path.dirname(currentRegList[0])
            for fileName in os.listdir(currentSubPath):
              if os.path.basename(fileName).find(searchPattern)!=-1:
                if os.path.basename(fileName).find(brightFeildString)!=-1:
                  brightFeildChannel1 = currentSubPath + '/' + fileName
                  brightFeildFound1   = 1
            if brightFeildFound1:
              self.downsampleAndRunAutoContrastInImageJ( brightFeildChannel1, 0 )
              subsampleFilename1 = os.path.splitext(brightFeildChannel1)[0]+'_subsample.tif'
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
                      if os.path.splitext(fileName)[1].find('tif')!=-1:
                        fullPathToFile = currentSubPath + '/' + fileName
                        round2FlourChannels.append( fullPathToFile )
                if brightFeildFound2:
                  #Start by downsampling images
                  self.downsampleAndRunAutoContrastInImageJ( brightFeildChannel2, 1 )
                  subsampleScaleFile = os.path.splitext(brightFeildChannel2)[0]+'.subsample'
                  f = open(subsampleScaleFile, 'r')
                  subsampleScaling = float(f.read())
                  f.close()
                  os.remove( subsampleScaleFile )
                  subsampleFilename2 = os.path.splitext(brightFeildChannel2)[0]+'_subsample.tif'
                  #Run gdbicp
                  args = [GdbIcpExec, subsampleFilename1, subsampleFilename2, '-model', '0', '-no_render', '1']
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
                            if line.find('AFFINE')!=-1 or counter1:
                              if line.find('AFFINE')!=-1:
                                f.close()
                              counter1 = counter1 + 1
                            if counter1==5:
                              stringArray = line.split()
                              xTranslate = -1*subsampleScaling*(float(stringArray[2])-float(stringArray[0]))
                              yTranslate = -1*subsampleScaling*(float(stringArray[3])-float(stringArray[1])) #Origin in ITK
                          if counter1<5:
                            print( 'Bad GDBICP file and skipping following file:' )
                            print currentRegList[i]
                          os.remove( xformFile )
                          TranslateExe = 'TranslateThenAffineRegister'+self.execExt
                          args = [ TranslateExe, brightFeildChannel1, brightFeildChannel2, str(xTranslate), str(yTranslate) ];
                          args.extend( round2FlourChannels );
                          self.RunExec( args, brightFeildChannel2 )
                          os.remove( brightFeildChannel2 )
##                          for removeName in round2FlourChannels:
##                            os.remove( removeName )
                else:
                  print( 'Failed to find the registration brightfield and skipping following file:' )
                  print currentRegList[i]
              os.remove( subsampleFilename1 )
            else:
              print( 'Failed to find the registration brightfield and skipping following files:' )
              print currentRegList
    print('I am done processing. Close the window.')

  #Run autocontras fiji/imagej macro
  def downsampleAndRunAutoContrastInImageJ( self, filename, writeScalingFlag ):
    ResampleExe = 'ResampleImage'+self.execExt
    if writeScalingFlag:
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
    nrrdConverter = self.execPref+'CziToNrrd'+self.execExt
    args = [ nrrdConverter, filename, self.numCoresToUse ]
    self.RunExec( args, filename )
    #Search and stitch new nrrd files into dicoms
    IlluminationEx = self.execPref+'Illumination'+self.execExt
    dicomConverter = self.execPref+'NrrdToDicom'+self.execExt
    filesThatExistPost = glob.glob(searchStr)
    for nrrdFile in filesThatExistPost:
      if nrrdFile.find(brightFieldStr0)!= -1:
        args = [ IlluminationEx, nrrdFile, '0', self.numCoresToUse, '1' ]
      else:
        if nrrdFile.find(Channel1Str) != -1 or nrrdFile.find(Channel3Str) != -1:
          args = [ IlluminationEx, nrrdFile, '38', self.numCoresToUse ] #Change 21
        else:
          args = [ IlluminationEx, nrrdFile, '0', self.numCoresToUse ]
      args1 = [ dicomConverter, nrrdFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
#      if nrrdFile.find(brightFieldStr0)!=-1:
#        self.RunExec( args1, nrrdFile ) #stitches raw uncorrected images
#      if nrrdFile.find(Channel3Str)!=-1 or nrrdFile.find(Channel1Str)!=-1:#nrrdFile.find(brightFieldStr0)!=-1 and 
      self.RunExec( args, nrrdFile )
      nrrdIllFile = os.path.splitext(nrrdFile)[0]+'_IlluminationCorrected.nrrd'
      if os.path.exists(nrrdIllFile):
        os.remove( nrrdFile )
        args = [ dicomConverter, nrrdIllFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
        self.RunExec( args, nrrdIllFile )
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
    '''f0 = open( errFilename, 'a' )
    f0.write( stderr )
    f0.close()'''

  def RunExec( self, args, filename ):
    '''errFilename = os.path.splitext(filename)[0]+'.errlog'
    logFilename = os.path.splitext(filename)[0]+'.log' '''
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    self.stdstrring += stdout
    self.stdsterr   += stderr
    '''print stdout, stderr
    f1 = open( errFilename, 'a' )
    f1.write( stderr )
    f1.close()
    f2 = open( logFilename, 'a' )
    f2.write( stdout )
    f2.close()'''

if __name__=='__main__':
#  PassFileNamesFromCLInsteadOfTk(sys.argv)  #No gui code
  root = Tkinter.Tk()
  TkFileDialog(root).pack()
  root.mainloop()

