import os, sys, glob, stat, getopt, platform, subprocess, multiprocessing, fnmatch
import Tkinter, Tkconstants, tkFileDialog

cziStr = '.czi'
cziFilterStr = '*.czi'
brightFieldStr0 = '_C0.nrrd'

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
    if numCores>14:
      numCores = 14
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
    print ('I am done. Close the window.')

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
      if nrrdFile.find(brightFieldStr0)!=-1:
        args = [ IlluminationEx, nrrdFile, self.numCoresToUse, '1' ]
      else:
        args = [ IlluminationEx, nrrdFile, self.numCoresToUse ]
      args1 = [ dicomConverter, nrrdFile, os.path.splitext(filename)[0]+'.xml', self.numCoresToUse ]
      self.RunExec( args1, nrrdFile ) #stitches raw uncorrected images
#      self.RunExec( args, nrrdFile )
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
    javaExec = 'C:\\Program Files\\Java\\jre7\\bin\\java.exe' #'java'+self.execExt
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
    errFilename = os.path.splitext(filename)[0]+'.errlog'
    logFilename = os.path.splitext(filename)[0]+'.log'
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

