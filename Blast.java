/* data structures */
import java.util.Hashtable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/* For file input and outputs */
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class Blast{
  //String db_filename;
  public static final int SLIDING_WINDOW_SIZE = 4;
  public Hashtable<String, String> db = new Hashtable<String, String>(); //when we first process the db.

  public String query;
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> q_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();
  //                kmer1  :{ entry1 : [1 , 3 , 6, 13], entry2 : [2, 7, 14, 21]}, kmer2... 
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> db_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();


  public double threshold = 1.0;

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

      //kmerizing query string.
      b.kmerize(b.q_kmer, SLIDING_WINDOW_SIZE, b.query, "query");

      System.out.println(b.q_kmer);

      for (String k : b.db.keySet()) {
        String val = b.db.get(k);
        b.kmerize(b.db_kmer, SLIDING_WINDOW_SIZE, val, k);
      }

      System.out.println(b.db_kmer);

      for (String key : b.q_kmer.keySet()) {
        b.matchAll(key, b.q_kmer.get(key).get("query"), b.db_kmer, b.threshold);  //..., b.q_kmer.get(key), ...
      }
      
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

  //ArrayList<ArrayList<String>>
  public void matchAll(String query_kmer, ArrayList<Integer> queryIndices, Hashtable<String, Hashtable<String, ArrayList<Integer>>> target, double threshold) {
    // for (String k : query.keySet()) {
    //   String val = query.get(k);
    // }

    Hashtable<String, ArrayList<Integer>> matches = target.get(query_kmer);

    // System.out.println("key: " + query_kmer); 
    
    if (matches != null) {
      // System.out.println("matches: " + matches.keySet());
      
      for (String dbName : matches.keySet()) {
        ArrayList<Integer> indices = matches.get(dbName);
        // System.out.println("match: " + dbName + " at: " + indices);
        //this.db

        batchExtendAlignment(query_kmer, queryIndices, dbName, this.db.get(dbName), indices);
      }
    }
    
  }

  public void batchExtendAlignment(String query_kmer, ArrayList<Integer> queryIndices, String db_key, String db, ArrayList<Integer> db_indices) {
    // System.out.println("query: " + this.query + "\nquerying kmer: " + query_kmer + "\nqueryIndices: " + queryIndices + "\ndb: " + db + "\ndb_indices: " + db_indices);
    for (Integer q_index : queryIndices) {
      for (Integer db_index : db_indices) {
        extendSingleAlignment(query_kmer, q_index.intValue(), db_key, db, db_index.intValue());
      }
    } 
  }

  public void extendSingleAlignment(String query_kmer, int q_index, String db_key, String db, int db_index) {
    System.out.println("in extendSingleAlignment");
    System.out.println("query: " + this.query + "\nquerying_kmer: " + query_kmer + "\nq_index: " + q_index + "\ndb: " + db + "\ndb_index: " + db_index);
    
    //backwards
    String back_qmatch = "";
    String back_dbmatch = ""; 
    int offset = 1;
    int j = db_index - offset;
    int back_i = q_index;
    for (int i = q_index - 1 ; i > -1; i--) {
      if (j > -1 && this.query.charAt(i) == db.charAt(j)) {
        back_qmatch = "" + this.query.charAt(i) + back_qmatch;
        back_dbmatch = "" + db.charAt(j) + back_dbmatch; 
        back_i = i;
      } else {
        break;
      }
      j = db_index-offset;
      offset++;
    }
    System.out.println("back-extended:");
    System.out.println(back_qmatch + " at: " + back_i + "\n" + back_dbmatch + " at: " + j);

    //forwards
    String forward_qmatch = "";
    String forward_dbmatch = "";
    
    offset = 0;
    j = db_index + query_kmer.length() + offset;

    int forward_i = q_index + query_kmer.length();
    for (int i = forward_i; i < this.query.length(); i++) {
      j = db_index + query_kmer.length() + offset;
      if (j < db.length() && this.query.charAt(i) == db.charAt(j)) {
        forward_qmatch += "" + this.query.charAt(i);
        forward_dbmatch += "" + db.charAt(j);
        forward_i++;
      } else {
        break;
      }
      offset++;
    }

    System.out.println("forward-extended:");
    System.out.println(forward_qmatch + " at: " + forward_i + "\n" + forward_dbmatch + " at: " + j);

    String extendedQuery = back_qmatch + query_kmer + forward_qmatch;
    String extendedDB = back_dbmatch + query_kmer + forward_dbmatch;

    System.out.println("results:");
    System.out.println(extendedQuery + "\n" + extendedDB);

  }

}