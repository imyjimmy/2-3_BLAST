/* data structures */
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.AbstractMap;

/* For file input and outputs */
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class Blast{
  //String db_filename;
  public Hashtable<String, String> db = new Hashtable<String, String>(); //when we first process the db.

  public String query;
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> q_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();
  //                kmer1  :{ entry1 : [1 , 3 , 6, 13], entry2 : [2, 7, 14, 21]}, kmer2... 
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> db_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();

  public void parseDB(String filename) {
    // System.out.println("public void parseFile(String filename) {");
    File file = new File(filename);
    try { 
        Scanner inputFile = new Scanner(file); 
        String key = "";
        String val = "";
        do {
          String line = inputFile.nextLine();
          if (line.startsWith(">")) {
            if (!"".equals(val) && !"".equals(key)) { //puts the previous key val pair into sequences.
              if (val.matches("[ATCG]+")) { // valid sequence. assuming no 'U' for now.
                db.put(key, val);
              } else { //invalid sequence
                System.out.println("This sequence contains terrible formatting \nand illegal sequences that came here illegally. \nMake your sequence great and try again.");
                System.exit(1);
              }
              key = "";
              val = "";
            } 
            key = line.substring(1);
          } else {
            val += line;
          }
        } while(inputFile.hasNextLine()); 
    } catch (FileNotFoundException fnfe) {
      fnfe.printStackTrace();
    }

    System.out.println(db);
  }

  public void parseQuery(String filename) {
    File file = new File(filename);
    try { 
        Scanner inputFile = new Scanner(file); 
        do {
          this.query = inputFile.nextLine();
          System.out.println("query: " + query);
          //...
        } while(inputFile.hasNextLine()); 
    } catch (FileNotFoundException fnfe) {
      fnfe.printStackTrace();
    }    
  }

  public static void main(String[] args) {
    Blast b = new Blast();
    Hashtable<String, String> parameters = new Hashtable<String, String>();
    
    if (args.length > 0 && (args[0].equals("-h") || args[0].equals("--help"))) {
      System.out.println("params: -h or --help to print this screen.");
      System.out.println("--database or -d [filepath]");
      System.out.println("--query or -q [filepath]");
      System.out.println("command line queries are not accepted, please put in a file.");
      System.exit(0);
    } else {
      for (int i = 0; i < args.length-1; i+=2) {
          parameters.put(args[i],args[i+1]);
      }    
      System.out.println(parameters);
      b.start(parameters);  

      //kmerizing.
      b.kmerize(b.q_kmer, 4, b.query, "query");

      System.out.println(b.q_kmer);

      for (String k : b.db.keySet()) {
        String val = b.db.get(k);
        b.kmerize(b.db_kmer, 4, val, k);
      }

      System.out.println(b.db_kmer);
    }
  }

  public void start(Hashtable<String, String> parameters) {
    String db_filename = parameters.get("-d");
    if (db_filename == null) {
      db_filename = parameters.get("--database");
      if (db_filename == null) {
        System.out.println("db file name is missing. Please supply -d [filepath] or --database [filepath]");
        System.exit(1);
      }
    }
    //this.db_filename = db_filename; //needed?
    String q_filename = parameters.get("-q");
    if (q_filename == null) {
      q_filename = parameters.get("--query");
      if (q_filename == null) {
        System.out.println("query file name is missing. Please supply -q [filepath] or --query [filepath]");
        System.exit(1);
      }
    }
    this.parseDB(db_filename);
    this.parseQuery(q_filename);
  }

  public void kmerize(Hashtable<String, Hashtable<String, ArrayList<Integer>>> kmer_db, int window, String toKmerize, String key) {
    for (int i = 0; i < toKmerize.length()-window+1; i++) {
      // System.out.println("substring to kmerize: " + toKmerize.substring(i, i+window));
      // kmer_db.put(toKmerize.substring(i, i+window), new Integer(i));
      String kmer = toKmerize.substring(i, i+window);
      Hashtable<String, ArrayList<Integer>> kmer_hits = kmer_db.get(kmer);
      ArrayList<Integer> indices;

      if (kmer_hits != null) {
        indices = kmer_hits.get(key);
        if (indices == null) {
          indices = new ArrayList<Integer>();
        }
        indices.add(new Integer(i));
      } else {
        kmer_hits = new Hashtable<String, ArrayList<Integer>>();
        indices = new ArrayList<Integer>();
        indices.add(new Integer(i));
      }

      kmer_hits.put(key, indices);
      kmer_db.put(kmer, kmer_hits);
    }
  
  }


  /*scratch */
        //   Map.Entry<String,Integer> pair1 = new AbstractMap.SimpleEntry<>("something",1);
        // Map.Entry<String,Integer> pair2 = new AbstractMap.SimpleEntry<>("something else",2);
        // b.db_kmerHits.add(pair1);
        // b.db_kmerHits.add(pair2);
        // // b.kmerize(b.db_kmer, 4);  
        // System.out.println(b.db_kmerHits);

}