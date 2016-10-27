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
  public static int sliding_window_size = 28;
  public double threshold = 1.0;

  public Hashtable<String, String> db = new Hashtable<String, String>(); //when we first process the db.

  public String query = "";
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> q_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();
  //                kmer1  :{ entry1 : [1 , 3 , 6, 13], entry2 : [2, 7, 14, 21]}, kmer2... 
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> db_kmer = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();

  // public Hashtable<String, Hashtable<Integer, ArrayList<Integer>>> db_hits = new Hashtable<String, Hashtable<Integer, ArrayList<Integer>>>();
  // db_hits
  // Hashtable<String, HashTable<Integer, String>>
  //           {"db1" => {3 => [12,21,45]}
  //           { db match => { db index => [a list of lengths]}}
  //           "db2": => {}
  //           }
  
  public Hashtable<String, Hashtable<String, ArrayList<Integer>>> hits = new Hashtable<String, Hashtable<String, ArrayList<Integer>>>();
  //              "db1" =>            {"3:12"   => [12,445,] }
  //              db match =>         { "db_index:length" => [query_indices] }
  
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

        // System.out.println("key: " + key + " val: " + val);
        db.put(key, val);
    } catch (FileNotFoundException fnfe) {
      fnfe.printStackTrace();
    }

    // System.out.println(db);
  }

  public void parseQuery(String filename) {
    File file = new File(filename);
    try { 
        Scanner inputFile = new Scanner(file); 
        while(inputFile.hasNextLine()) {
          String line = inputFile.nextLine();
          if (line != null) {
            this.query += line;
          }
          //...
        }
        System.out.println("query: " + query);
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
      System.out.println("--threshold or -t [1.0 or less]");
      System.out.println("--kmer or -k [number]");
      System.out.println("command line queries are not accepted, please put in a file.");
      System.exit(0);
    } else {
      for (int i = 0; i < args.length-1; i+=2) {
          parameters.put(args[i],args[i+1]);
      }    
      System.out.println(parameters);
      b.start(parameters);  

      //kmerizing query string.
      b.kmerize(b.q_kmer, sliding_window_size, b.query, "query");

      // System.out.println(b.q_kmer);

      for (String k : b.db.keySet()) {
        String val = b.db.get(k);
        b.kmerize(b.db_kmer, sliding_window_size, val, k);
      }

      // System.out.println(b.db_kmer);

      for (String key : b.q_kmer.keySet()) {
        b.matchAll(key, b.q_kmer.get(key).get("query"), b.db_kmer, b.threshold);  //..., b.q_kmer.get(key), ...
      }
      
      System.out.println(b.hits);
      b.mergeHits();
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

    String threshold = parameters.get("-t");
    if (threshold == null) {
      threshold = parameters.get("--threshold");
      if (threshold == null) {
        System.out.println("threshold is missing, using the default 1.0");
        threshold = "1.0";
        // System.exit(1);
      }
    }

    this.threshold = Double.parseDouble(threshold);

    String kmer = parameters.get("-k");
    if (kmer == null) {
      kmer = parameters.get("--kmer");
      if (kmer == null) {
        System.out.println("kmer is missing, defaulting to 28");
        kmer = "28";
      }
    }
    this.sliding_window_size = Integer.parseInt(kmer);

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

  //return ArrayList<ArrayList<String>> (?)
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
    // System.out.println("in extendSingleAlignment");    
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
      offset++;
      j = db_index-offset;
    }
    // System.out.println("back-extended:");
    // System.out.println(back_qmatch + " at: " + back_i + "\n" + back_dbmatch + " at: " + j);
    int back_q_index = back_i;
    int back_db_index = j+1;

    //forwards
    String forward_qmatch = "";
    String forward_dbmatch = "";
    
    offset = 1;
    j = db_index + query_kmer.length() + offset;

    int forward_i = q_index + query_kmer.length();
    for (int i = forward_i; i < this.query.length(); i++) {
      j = db_index + query_kmer.length() + offset - 1;
      if (j < db.length() && this.query.charAt(i) == db.charAt(j)) {
        forward_qmatch += "" + this.query.charAt(i);
        forward_dbmatch += "" + db.charAt(j);
        forward_i++;
      } else {
        break;
      }
      offset++;
    }

    // System.out.println("forward-extended:");
    // System.out.println(forward_qmatch + " at: " + forward_i + "\n" + forward_dbmatch + " at: " + j);

    String extendedQuery = back_qmatch + query_kmer + forward_qmatch;
    String extendedDB = back_dbmatch + query_kmer + forward_dbmatch;

    //todo: extend query past exact matches, up to threshold.
    System.out.println("querying_kmer: " + query_kmer + "\nq_index: " + back_q_index + "\ndbname: " + db_key + "\ndb_index: " + back_db_index);
    System.out.println("results:");
    System.out.println(extendedQuery + "\n" + extendedDB);
    if (! (this.threshold == 1.0)) {
      System.out.println("extending query to threshold.");
      extendToThreshold(extendedQuery, back_q_index, db_key, back_db_index);
    }

  }

  public void extendToThreshold(String match, int q_index, String db_key, int db_index) {
    String back_qmatch = "";
    String back_dbmatch = "";
    boolean back_end = false;

    String forward_qmatch = "";
    String forward_dbmatch = "";
    boolean forward_end = false;

    double match_fraction = 1.0;
    int i = 0;
    int num_mismatches = 0;
    int front_i = 0;
    int back_i = 0;

    String db = this.db.get(db_key);


    while (match_fraction > this.threshold) {
      // System.out.println("inside while");

      if (i % 2 == 0 && !back_end) { //back way
        // System.out.println("back way");
        //q_index - 1 - i  vs db_index - 1 - i
        if (q_index-1-back_i > -1 && db_index-1-back_i > -1) {
          if (this.query.charAt(q_index-1-back_i) == db.charAt(db_index-1-back_i)) {
            //nothing to actually do here
          } else {
            num_mismatches++;
            match_fraction = 1.0 - ((double) num_mismatches / (double) (match.length()+forward_qmatch.length()+back_i+1));
          }

          back_qmatch = "" + this.query.charAt(q_index-1-back_i) + back_qmatch;
          back_dbmatch = "" + db.charAt(db_index-1-back_i) + back_dbmatch;
          back_i++;
        } else {
          //do not pass go, do not collect 200
          // System.out.println("reached the back end. (start of string)");
          back_end = true;
        }

      } else if (i % 2 == 1 && !forward_end) { //front way
        //q_index+match.length() + i vs db_index + match.length() + i
        // System.out.println("going front way");
        if (q_index+match.length()+front_i < this.query.length() && db_index + match.length() + front_i < db.length()) {
          int q = q_index+match.length()+front_i;
          int d = db_index+match.length()+front_i;
          if (this.query.charAt(q_index+match.length()+front_i) == db.charAt(db_index+match.length()+front_i)) {
            //
            // System.out.println("there was a match");
          } else {
            // System.out.println("mismatch");
            // System.out.println("at query: " + q + " at db: " + d + " comparing: " + this.query.charAt(q_index+match.length()+front_i) + " with " + db.charAt(db_index+match.length()+front_i));
            num_mismatches++;
            match_fraction = 1.0 - ((double) num_mismatches / (double) (match.length() + back_qmatch.length() + front_i + 1));
            // System.out.println("match_fraction: " + match_fraction);
          }

          forward_qmatch += "" + this.query.charAt(q_index+match.length()+front_i);
          forward_dbmatch += "" + db.charAt(db_index+match.length()+front_i);
          front_i++;

          // System.out.println("query forward: " + forward_qmatch + "\ndb forward: " + forward_dbmatch);
        } else {
          //do not pass go...
          // System.out.println("reached the forward end (end of string)");
          forward_end = true;
        }

      } else {
        if (back_end && forward_end) {
          // System.out.println("breaking...");
          break;
        }
        // System.out.println("wierd else place");
      }

      i++;
    }

    q_index = q_index - back_qmatch.length();
    db_index = db_index - back_dbmatch.length();

    String extendedQuery = back_qmatch + match + forward_qmatch;
    String extendedDB = back_dbmatch + match + forward_dbmatch;

    // System.out.println("q_index: " + q_index + "\ndbname: " + db_key + "\ndb_index: " + db_index);

    // System.out.println("results:");
    // System.out.println(extendedQuery + "\n" + extendedDB);

    addHit(extendedQuery, q_index, db_key, extendedDB, db_index);
    printConsensus(extendedQuery, q_index, db_key, extendedDB, db_index);
  }

  public void addHit(String query, int q_index, String db_key, String db_substr, int db_index) {
    //hits:   
    //              "12:21" =>            {"db1"   => [12,445,] }
    //              beg:end =>         { "db_name" => [query_indices] }

    int endpos = q_index + query.length() - 1;
    String q_key = Integer.toString(q_index) + ":" + Integer.toString(endpos);
    // System.out.println(q_key); 

    Hashtable<String, ArrayList<Integer>> hits_q_key = hits.get(q_key);
    if (hits_q_key == null) {
      hits_q_key = new Hashtable<String, ArrayList<Integer>>();
    }

    ArrayList<Integer> indices = hits_q_key.get(db_key);
    if (indices == null) {
      indices = new ArrayList<Integer>();
    }

    if (!indices.contains(db_index)) {
      indices.add(db_index);
    } else {
      System.out.println("already contains this index.");
    }

    hits_q_key.put(db_key, indices);
    hits.put(q_key, hits_q_key);

  }

  public void mergeHits() {
    ArrayList<String> keys = new ArrayList<String>(hits.keySet());
    for (int i = 0; i < keys.size(); i++) {
      String key = keys.get(i);
      String[] beg_end = key.split(":");
      int begin = Integer.parseInt(beg_end[0]);
      int end = Integer.parseInt(beg_end[1]);

      ArrayList<String> keys_in_range = getKeysInRange(begin, end, key);

      Hashtable<String, ArrayList<Integer>> hits_of_coord = hits.get(key);
      ArrayList<String> dbs = new ArrayList<String>(hits_of_coord.keySet());

      for (String key_in_r : keys_in_range) {
        String[] r_beg_end = key_in_r.split(":");
        int targetDistance = begin - Integer.parseInt(r_beg_end[0]); //looking for this distance.
        //          -4     =    12 - 16
        Hashtable<String, ArrayList<Integer>> hits_in_range = hits.get(key_in_r);
        for (String db : dbs) { //a particular db name that appears in dbs, the hashtable for a given query key.
          ArrayList<Integer> hits_in_merge_range = hits_in_range.get(db);
          ArrayList<Integer> hits_in_db = hits_of_coord.get(db);
          if (hits_in_merge_range == null) {
            //remove...or skip...
            //keys_in_range.remove(key_in_r);
          } else { //hits in db is a thing.
            //loook for a target index a certain distance away...
            //begin - targetDistance
            for (Integer _i : hits_in_db) {
              int targetIndex = _i.intValue() - targetDistance;
              if (hits_in_merge_range.contains(targetIndex)) {
                //merge.
                //meaning...? key, db, _i, key_in_r, _i.intValue() - target
              }
            }
            
          }          
        }

      }
      
      // for (String db : dbs) {
      //   System.out.println(db);
      //   ArrayList<Integer> match_indices = hits_of_coord.get(db);
      // }
    }
  }

  public ArrayList<String> getKeysInRange(int begin, int end, String exclude) {
    ArrayList<String> keys = new ArrayList<String>(hits.keySet());

    ArrayList<String> toReturn = new ArrayList<String>();
    for (String key : keys) {
      if (key.compareTo(exclude) != 0) {
        String[] beg_end = key.split(":");
        int _begin = Integer.parseInt(beg_end[0]);
        int _end = Integer.parseInt(beg_end[1]);
        if ( begin <= _begin && end >= _begin) { //sequence to merge occurs after
          toReturn.add(key);
        } else if ( _begin <= begin && _end >= begin) { //before
          toReturn.add(key);
        }
      } else {
        //skip
      }
    }

    return toReturn;
  }

  public void printConsensus(String query, int q_index, String db_key, String db_substr, int db_index) {
    String buffer = "";
    int db_index_parse = db_index;
    while (db_index_parse / 10 != 0) {
      buffer += " ";
      db_index_parse = db_index_parse / 10;
    }
    if (query.length() < 60) {
      System.out.print(buffer + q_index);
      for (int i = 0; i < query.length()-2; i++) {
        System.out.print(" ");
      }
      int q_end = q_index + query.length()-1;
      System.out.println(q_end + " query");

      //query string
      System.out.println(buffer + query);

      //horiz hashes
      System.out.print(buffer);
      for (int i = 0; i < query.length(); i++) {
        System.out.print("|");
      }
      System.out.println();
      
      //db string
      System.out.println(buffer+db_substr);

      //db index
      System.out.print(db_index);
      for (int i = 0; i < db_substr.length()-2; i++) {
        System.out.print(" ");
      }
      int db_end = db_index + db_substr.length() - 1;
      System.out.println(db_end + " " + db_key);

      //----//
      System.out.print(buffer);
      for (int i = 0; i < query.length(); i++) {
        System.out.print("-");
      }
      System.out.println();

      //concensus
      System.out.print(buffer);
      for (int i = 0; i < query.length(); i++) {
        if (query.charAt(i) == db_substr.charAt(i)) {
          System.out.print("*");
        } else {
          System.out.print("$");
        }
      }

      System.out.println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    } else {

      //query, q_index, db_key, db_substr
      System.out.println("query match at: " + q_index + " with db: " + db_key + " at: " + db_index); 
      int bottom_index = 0;
      for (int i = 0; i < query.length(); i++) {             
          System.out.print(query.charAt(i));    
          if (((i % 60 == 0) && i > 0) || (i == query.length()-1)) {
              System.out.print("\n");
              for (int j = bottom_index; j<=i; j++) {
                  System.out.print("|");
                  if (j == i) {
                      System.out.print("\n");
                  }
              }
              for (int k = bottom_index; k<=i; k++) {
                  System.out.print(db_substr.charAt(k));  
              }
              System.out.println();
              for (int l = bottom_index; l<=i; l++) {
                System.out.print("-");
              }
              System.out.println();
              for (int m = bottom_index; m<=i; m++) {
                if (db_substr.charAt(m) == query.charAt(m)) {
                  System.out.print("*");
                } else {
                  System.out.print("$");
                }
                bottom_index++;
              }
              System.out.println("\n");
          }
      } 
      System.out.print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }

  }




}