/*
 * 
 * Wave wrote the script using the classes from picard... 
 * picard and sam 1.48 version
 * 
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.SolexaQualityConverter;
import net.sf.samtools.util.StringUtil;

/**
 * De-barcode the sequences into different samples...
 * See <a href="http://maq.sourceforge.net/fastq.shtml">MAQ FastQ specification</a> for details.
 * Three fastq versions are supported: FastqSanger, FastqSolexa and FastqIllumina.
 * Input files can be in GZip format (end in .gz).
 */
public class deBCSCell extends CommandLineProgram {
    private static final Log LOG = Log.getInstance(deBCSCell.class);

    @Usage 
    public String USAGE = "\tdeBCSCell V1.0.\n" +
    		"\tExtract the barcod sequences and sign the correspoding sequences to samples....\n"
        + "\tInput files can be in GZip format (end in .gz).\n\n" +
        		"Usage: \n"
//        +"\tjavasol deBarcoding F1=seq1.txt.gz BL=sampleBarcode.txt \n" 
        		+"\tjavasol deBarcoding F1=seq1.txt.gz F2=bcseq2.txt.gz BL=sampleBarcode.txt"
        ;


    @Option(shortName="F1", doc="Input fastq file (optionally gzipped) for single end data")
    public File FASTQ;

    @Option(shortName="F2", doc="Input fastq file (optionally gzipped) for barcode sequences.", optional=true)
    public File FASTQ2;

    @Option(shortName="BL", doc="Barcode sequence list and sample name...")
    public File bcFile;

    @Option(shortName="DR", doc="Directory for the sequences...", optional=true)
    public File direc=new File("sequences/");

    @Option(shortName="ST", doc="Options to do stat ONLY... \nFalse (default) for stat and write out sequences;\ntrue for stat only ...\n", optional=true)
    public boolean stat=false;
    
    @Option(shortName="NS", doc="Number of sequences for different barcodes...", optional=true)
    public File numStat=new File("numStat.txt");

    @Option(shortName = "MD", doc = "Max mismatch to the sample barcodes..", optional = true)
    public static int maxDis=1; 
    
    public boolean randomBarcodes=false;
    public int ranBarLen=6;
    public int samBarLen=8;
    

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    /** Stock main method. */
    public static void main(final String[] argv) {
    	System.out.print(argv[1]);
        System.exit(new deBCSCell().instanceMain(argv));
        
    }
    
    /* Simply invokes the right method for unpaired or paired data. */
    protected int doWork() {
        final int readCount = (FASTQ2 == null) ?  doUnpaired() : doPaired();
        LOG.info("Processed " + readCount + " fastq reads");
//        String[] xx=allBarSeq();
//        System.out.println(xx);
        return 0;
    }

   
    protected int doUnpaired() {
    	
        IoUtil.assertFileIsReadable(FASTQ);
//        IoUtil.assertFileIsWritable(OUTPUT);
        
        final FastqReader freader = new FastqReader(FASTQ);
        System.out.print("SHOULD BE TWO FASTQ FILES.\n");
        System.out.print("Now is running the remove molecular barcodes and mult-G's in the beginnning of reads...\n");
        
        int readCount=0;
        String newFile= FASTQ.getName().replace(".fastq.gz", ".fastq");
        FastqWriter wr = new FastqWriter(new File(newFile));
        
        for ( ; freader.hasNext()  ; readCount++) {
            final FastqRecord frec1 = freader.next();
           
          String barSeq1=frec1.getReadString().substring(0,ranBarLen);

          String seq1=frec1.getReadString().substring(ranBarLen+1);
//          System.out.println(frec1.getReadHeader());

          int ii=0;
          for (ii = 0; (seq1.charAt(ii)=='G' | seq1.charAt(ii)=='g' | seq1.charAt(ii)=='N' | seq1.charAt(ii)=='n') & ii < seq1.length()-1  ; ii++){}
          
          if(ii<14){
        	  
        	  String seqn=frec1.getReadString().substring(ranBarLen+ii+1);
              String quan=frec1.getBaseQualityString().substring(ranBarLen+ii+1);
//              System.out.println(frec1.getReadHeader());
              String readnameN=frec1.getReadHeader().replaceAll(" ", "_")+":"+barSeq1+":"+ii;
              
              FastqRecord nfrec1= new FastqRecord(readnameN, seqn,
            		  readnameN, quan);
              
             
              wr.write(nfrec1);
               
          }else{
        	  System.out.println(seq1+":"+ii);
          }
           
          if( (readCount % 1000000 ) ==0 & (readCount / 1000000) >0){
          	System.out.println(readCount / 1000000 + " Million reads are processed");
          	}
          
          
        }
        wr.close();
        return 0;
    }

    /** More complicated method that takes two fastq files. */
    protected int doPaired() {
        IoUtil.assertFileIsReadable(FASTQ);
        IoUtil.assertFileIsReadable(FASTQ2);
//        IoUtil.assertFileIsWritable(OUTPUT);
        
        HashMap wl=readBarcodeList(stat);
//        
        // stupid, should change...
        LinkedHashMap num=new LinkedHashMap();
        Iterator it =wl.keySet().iterator();
        while(it.hasNext()){
    		String ii=it.next().toString();
    		num.put(ii, 0);
    	}
    	
//    	
    	//
        HashMap mpbarcode=getBarcode(wl,samBarLen);
        
        System.out.println("Debarcoding is processing...");
        
        final FastqReader freader1 = new FastqReader(FASTQ);
        final FastqReader freader2 = new FastqReader(FASTQ2);
        
        int readCount = 0;
        
        for ( ; freader1.hasNext() && freader2.hasNext() ; readCount++) {
            final FastqRecord frec1 = freader1.next();
            final FastqRecord frec2 = freader2.next();
           
          String barSeq1=frec1.getReadString().substring(0,ranBarLen);

          String seq1=frec1.getReadString().substring(ranBarLen+1);
//          System.out.println(frec1.getReadHeader());

          int ii=0;
          for (ii = 0; (seq1.charAt(ii)=='G' | seq1.charAt(ii)=='g' | seq1.charAt(ii)=='N' | seq1.charAt(ii)=='n') & ii < seq1.length()  ; ii++){}
          
          
          String seqn=frec1.getReadString().substring(ranBarLen+ii+1);
          String quan=frec1.getBaseQualityString().substring(ranBarLen+ii+1);
//          System.out.println(frec1.getReadHeader());
          String readnameN=frec1.getReadHeader().replaceAll(" ", "_")+":"+barSeq1+":"+ii;
          
          FastqRecord nfrec1= new FastqRecord(readnameN, seqn,
        		  readnameN, quan);
          
          String barSeq2=frec2.getReadString();
//        System.out.println(barSeq);
          
          
          if(mpbarcode.containsKey(barSeq2)){
          	String idBarcode=mpbarcode.get(barSeq2).toString();
              
          	if(!stat){
          		
//                System.out.println(barSeq+"\t"+idBarcode+"\t"+(barSeq.equalsIgnoreCase(idBarcode)));
                
                FastqWriter wr= (FastqWriter) wl.get(idBarcode);
               
                if(wr!=null){
                	wr.write(nfrec1);
                	
                }
          	}

//            Stupid ...
            	int barcodeNum = Integer.parseInt(num.get(idBarcode).toString()) + 1;
            	num.put(idBarcode, barcodeNum);
//            
          	
             
          }else{
//          	System.out.println(barSeq+"\t");
          }
 
          
          
          /*
         * Print out the nubmer of reads processed...
         * */
        if( (readCount % 1000000 ) ==0 & (readCount / 1000000) >0){
        	System.out.println(readCount / 1000000 + " Million reads are processed");
        	}
        
        }
        



        if (freader1.hasNext() || freader2.hasNext()) {
            throw new PicardException("Input paired fastq files must be the same length");
        }
        closeWriter(wl);
        
        
        try {
			PrintWriter out1=new PrintWriter(new BufferedWriter(new FileWriter(numStat)));
			Iterator itt=num.keySet().iterator();
			while(itt.hasNext()){
	    		String ii=itt.next().toString();
//	    		System.out.println(ii+"\t"+ lhm.get(ii));
	    		out1.println(ii+"\t"+ num.get(ii));
	    	}
			out1.close();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        return readCount;
    }


    /** Returns read baseName and asserts correct pair read name format:
     * <ul>
     * <li> Paired reads must either have the exact same read names or they must contain at least one "/"
     * <li> and the First pair read name must end with "/1" and second pair read name ends with "/2"
     * <li> The baseName (read name part before the /) must be the same for both read names
     * <li> If the read names are exactly the same but end in "/2" or "/1" then an exception will be thrown 
     * </ul>
     */
    String getBaseName(final String readName1, final String readName2, final FastqReader freader1, final FastqReader freader2) {
        String [] toks = getReadNameTokens(readName1, 1, freader1);
        final String baseName1 = toks[0] ;
        final String num1 = toks[1] ;

        toks = getReadNameTokens(readName2, 2, freader2);
        final String baseName2 = toks[0] ;
        final String num2 = toks[1];

        if (!baseName1.equals(baseName2)) {
            throw new PicardException(String.format("In paired mode, read name 1 (%s) does not match read name 2 (%s)", baseName1,baseName2));
        }

        final boolean num1Blank = StringUtil.isBlank(num1);
        final boolean num2Blank = StringUtil.isBlank(num2);
        if (num1Blank || num2Blank) {
            if(!num1Blank) throw new PicardException(error(freader1,"Pair 1 number is missing (" +readName1+ "). Both pair numbers must be present or neither."));       //num1 != blank and num2   == blank
            else if(!num2Blank) throw new PicardException(error(freader2, "Pair 2 number is missing (" +readName2+ "). Both pair numbers must be present or neither.")); //num1 == blank and num =2 != blank 
        } else {
            if (!num1.equals("1")) throw new PicardException(error(freader1,"Pair 1 number must be 1 ("+readName1+")"));
            if (!num2.equals("2")) throw new PicardException(error(freader2,"Pair 2 number must be 2 ("+readName2+")"));
        }

        return baseName1 ;
    }

    /** Breaks up read name into baseName and number separated by the last / */
    private String [] getReadNameTokens(final String readName, final int pairNum, final FastqReader freader) {
        if(readName.equals("")) throw new PicardException(error(freader,"Pair read name "+pairNum+" cannot be empty: "+readName));

        final int idx = readName.lastIndexOf("/");
        final String result[] = new String[2];

        if (idx == -1) {
            result[0] = readName;
            result[1] = null;
        } else {
            result[1] = readName.substring(idx+1, readName.length()); // should be a 1 or 2
            
            if(!result[1].equals("1") && !result[1].equals("2")) {    //if not a 1 or 2 then names must be identical
                result[0] = readName;
                result[1] = null;
            }
            else {
                result[0] = readName.substring(0,idx); // baseName
            }
        }

        return result ;
    }
    private LinkedHashMap readBarcodeList(boolean stat){
    	LinkedHashMap wm=new LinkedHashMap();
    	BufferedReader br = new BufferedReader(new InputStreamReader(IoUtil.openFileForReading(bcFile)));
    	try {
			while(br.ready()){
				String str=br.readLine();
				String a[]=str.split("\t");
				String sample=a[0];
				String bseq=a[1];
				
				
		    	if(bseq.length()!=samBarLen){
					System.out.println(bseq+":"+bseq.length()+"\tWarning: Barcode length is not equal to BARLEN: "+ samBarLen + "bp!! The first "+Math.min(bseq.length(), samBarLen)+" sequences will be compared....");
				}
		    	
				//check the directory
				if(direc.isDirectory() | stat){
//					System.out.println(direc);
					}else{
//						System.out.println("No directory called "+direc);
						direc.mkdirs();
				}
				
//				if(FASTQ2 == null){ //single end
					
					String fileName=direc+"/"+sample+"_sequence.txt";
					
					if(stat){
						wm.put(bseq, bseq);
					}else{
						
						FastqWriter writer = new FastqWriter(new File(fileName));
					
						wm.put(bseq,writer);
					}
					
//				}else{//paired end
//
//					String fileName1=direc+"/"+sample+"_p_1_sequence.txt";
//					String fileName2=direc+"/"+sample+"_p_2_sequence.txt";
//					
//					if(stat){
//						
//						wm.put(bseq, bseq);
//					
//					}else{
//						FastqWriter writer1 = new FastqWriter(new File(fileName1));
//						FastqWriter writer2 = new FastqWriter(new File(fileName2));
//						LinkedHashMap hmm=new LinkedHashMap();
////						System.out.println(hmm);
//						hmm.put(1,writer1);
//						hmm.put(2,writer2);
//						
//						wm.put(bseq,hmm);
//					}
//					}
				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print(e);
		}
    	return(wm);
    }
    
    private void closeWriter(HashMap hm){
//    	System.out.println(hm.keySet());
    	List list = new ArrayList(hm.keySet());
    	Iterator ite =list.iterator();
    	while(ite.hasNext()){
    		String ii=ite.next().toString();
    		
    		if(stat){}else{
    			
//    			if(FASTQ2 == null){
        			FastqWriter writer=(FastqWriter) hm.get(ii);
//            		System.out.println(ii);
            		writer.close();
//        		}else{
//        			FastqWriter writer1=(FastqWriter) ((HashMap) hm.get(ii)).get(1);
//        			FastqWriter writer2=(FastqWriter) ((HashMap) hm.get(ii)).get(2);
////            		System.out.println(ii);
//            		writer1.close();
//            		writer2.close();
//        		}
        		
    		}
    		
    	}
    	
    }
    
    
    /*
     * Creating all the possible barcodes in given length (barcodeLen)...
     *  5^barcodeLen
     * */
    private static String[] allBarSeq(int barcodeLen){
    	String[] abs=new String[(int)Math.pow(5,barcodeLen)];
    	
    	String[] bs=new String[]{"A","T","G","C","N"};
    	
    int[] byt=new int[barcodeLen];
//    System.out.println(barcodeLen);
    
    for(int i=0;i<abs.length;i++){
    		String sn= Integer.toString( i,  5);
//    		System.out.println(sn);
    		
    		for(int ii=0;ii<sn.length();ii++){
//    			System.out.println(barcodeLen-(sn.length()-i));
    			byt[barcodeLen-(sn.length()-ii)]=Integer.parseInt(sn.substring(ii,ii+1));
    		}
    		
    		StringBuffer seq=new StringBuffer();
    		for(int j=0;j<barcodeLen;j++){
    			seq=seq.append(bs[byt[j]]);
    		}
    		
    		abs[i]=seq.toString();
//    		System.out.println(abs[i]);
    		
    }    	
    	
//    	System.out.println(abs.length);
    System.out.println("Array for "+barcodeLen+" bp barcode sequences, with "+abs.length+" possibilities was created....");
    
    	return(abs);
    }

    private static double distance(String que, String sbj){
    	double dis=0;
    	
//    	if(que.length()!=sbj.length()){
//			System.out.println("Warning: Barcode length is not equal to BARLEN: "+ sbj.length() + "bp!! The first "+Math.min(que.length(), sbj.length())+" sequences will be compared....");
//		}
    		
    	for(int i=0;i<que.length();i++){
    		double dd=0;
    		if(que.substring(i,i+1).equalsIgnoreCase( sbj.substring(i,i+1))){
    			dd=0;
    		}else{
    			if(que.substring(i,i+1).equalsIgnoreCase("N")|sbj.substring(i,i+1).equalsIgnoreCase("N")){
    				dd=0.75;
    				}else{dd=1;}
    		}
    		dis=dis+dd;
    	}
    	return (dis);
    }
  
    private static HashMap getBarcode(HashMap vm,int barcodeLen){
    	
    	HashMap hm=new HashMap();
    	
    	String[] abs=allBarSeq(barcodeLen);
    	
    	Set key=vm.keySet();
    	String[] bseq =(String[]) key.toArray(new String[key.size()]);
    	
    	System.out.println("Assigning all possible barcode sequences to real barcode list....");

    	double[][] ad=new double[abs.length][bseq.length];
    	for(int i=0; i<abs.length;i++){
    		double min=0;
    		for(int j=0;j<bseq.length;j++){
    			
    			String que=bseq[j];
    			String sbj=abs[i];
    			double dd=distance(que,sbj);
    			ad[i][j]=dd;
    			
    			if(j==0){min=dd;}else{min=Math.min(min,dd);}
    		}
    		
    		int jj=0;
    		int jnum=0;
    		if(min<=maxDis){
    			for(int j=0;j<vm.size();j++){
    				if(ad[i][j]==min){
    					jj=j;jnum=jnum+1;
    				}
    			}
    			if(jnum==1){
    				hm.put(abs[i], bseq[jj]);
//    				System.out.println(abs[i]+"="+bseq[jj]);
    			}
    		}
    		
    		
    	}
    	
    	return(hm);
    }
    
    /** Little utility to give error messages corresponding to line numbers in the input files. */
    private String error(final FastqReader freader, final String str) {
        return str +" at line "+freader.getLineNumber() +" in file "+freader.getFile().getAbsolutePath();
    }

    // Read names cannot contain blanks
    private String getReadName(final String fastaqHeader) {
        final int idx = fastaqHeader.indexOf(" ");
        return (idx == -1) ? fastaqHeader : fastaqHeader.substring(0,idx); 
    }
}

