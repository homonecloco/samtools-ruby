$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
$: << File.expand_path('.')
require "test/unit"
require "bio/db/sam"
require "bio/db/sam/sam"


class TestBioDbSam < Test::Unit::TestCase

  #Set up the paths
  def setup
    @test_folder                = "test/samples/small"
    @testTAMFile                = @test_folder + "/test.tam"
    @testBAMFile                = @test_folder + "/testu.bam"
    @testReference              = @test_folder + "/test_chr.fasta"
    
  end
  
  #Removing the index files
  def teardown
    begin
      File.delete(@testReference + ".fai")
      p "deleted: " + @testReference + ".fai "
    rescue
    end
    begin
      File.delete(@testBAMFile + ".fai")
      p "deleted: " + @testBAMFile + ".bai "
    rescue
    end
  end


  #Test just opens and close a sam file for reading, from the low level API. 
  def test_openSAMFile
    bamfile                      = Bio::DB::SAM::Tools.samopen(@testTAMFile,"r",nil)
    Bio::DB::SAM::Tools.samclose(bamfile)
    assert(true, "file open and closed")
  end

  #Tests that an Exception is thrown when no parameters are given. 
  def test_new_class_empty
    begin
      bam                        = Bio::DB::Sam.new({})
      assert(false, "Should fail while opening without parameters")
    rescue Bio::DB::SAMException => e
      puts e.message
      assert(true, e.message)
    end
  end

  #Tests that an Exception is thrown when the bam file doesn't exists
  def test_new_class_empty_invalid_path
    begin
      sam                        = Bio::DB::Sam.new({:bam=>"INVALID"})
      sam.open
      sam.close
      assert(false, "Should fail with an invalid path")
    rescue Bio::DB::SAMException => e
      puts e.message
      assert(true, e.message)
    end
  end

  #Testst that you can open and close without a fasta file. 
  def test_class_text_read_no_faidx
    sam                          = Bio::DB::Sam.new({:tam=>@testTAMFile})
    sam.open
    sam.close
    assert(true, "file open and closed with the class")
  end
  
  #Tests what happens when a bam file is garbage collected. The object is destroyed. Used with flags
  #it shows that the file is closed prior the garbage colletion. 
  def test_class_text_read_no_close

    fam                          = Bio::DB::Sam.new({:tam=>@testTAMFile})
    fam.open
    fam                          = nil   
    ObjectSpace.garbage_collect

    assert(true, "file openend but not closed")
  end

  #Tests what happens when a tam file is garbage collected. The object is destroyed. Used with flags
   #it shows that the file is closed prior the garbage colletion.
  def test_class_binary_read_no_close
    Bio::DB::Sam.new({:bam=>@testBAMFile}).open
    ObjectSpace.garbage_collect
    assert(true, "BINARY file openend but not closed")  
  end

  #Tests that the coverage is readed for the whole chromosme. 
  def test_read_coverage
     sam       = Bio::DB::Sam.new({:bam=>@testBAMFile, :fasta=>@testReference})
     sam.open
  	File.open( @test_folder +"/ids2.txt", "r") do |file|
  	  puts "file opened"
  	  file.each_line{|line|
        fetching = line.split(' ')[0]
        puts "fetching: " + fetching
  	    sam.load_reference
    	  seq = sam.fetch_reference(fetching, 0, 16000)
  	    als = sam.fetch(fetching, 0, seq.length) 
     #      p als
  	  if als.length() > 0 then
  	       p fetching
                p als
  	    end
      }

  	end
    sam.close
   assert(true, "Finish")
  end
  

#  def test_read_TAM_as_BAM
#    begin
#      sam                          = Bio::DB::Sam.new({:bam=>@testTAMFile})
#      sam.open
#      sam.close
#      assert(false, "Should raise an exception for reading a BAM as TAM") 
#    rescue Bio::DB::SAMException => e
#      assert(true, "Properly handled")
#    end 
#  end

# def test_read_BAM_as_TAM
#    begin
#      sam                          = Bio::DB::Sam.new({:tam=>@testBAMFile})
#      sam.open
#      sam.close
#      assert(false, "Should raise an exception for reading a BAM as TAM") 
#    rescue Bio::DB::SAMException => e
#      assert(true, "Properly handled")
#    end 
#  end

#Tests that the index is created
  def test_bam_load_index
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    index = sam.load_index
    sam.close
    assert(true, "BAM index loaded")
  end

#Tests that the TAM file can't be used, since it can't be indexed
  def test_tam_load_index
    begin
      sam       = Bio::DB::Sam.new({:tam=>@testTAMFile})
      sam.open
      sam.load_index
      sam.close
      assert(false, "TAM index loaded")
    rescue Bio::DB::SAMException => e
      assert(true, "Unable to load an index for a TAM file")
    end
  end

#Tests that a segment fo the chromosme is readed correctly. 
  def test_read_segment
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    als = sam.fetch("chr_1", 0, 500)
    p als 
    sam.close
    assert(true, "Seems it ran the query")
    #node_7263       238     60 has 550+, query from 0 to 500, something shall come.... 
  end
  
  
  #Testst that you can fetch a segment from a preparsed index. 
  def test_read_segment_from_index
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    ind = sam.get_chr_index("chr_1")
    als = sam.fetch(ind, 0, 500)
    p als 
    sam.close
    assert(true, "Seems it ran the query")
    #node_7263       238     60 has 550+, query from 0 to 500, something shall come.... 
  end

  #Tests that an exception is thrown when the chromosome is not a valid chromosome. 
  def test_read_invalid_reference
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    begin
      als = sam.fetch("Chr1", 0, 500)
      p als 
      sam.close
      assert(false, "Seems it ran the query")
    rescue Bio::DB::SAMException => e
      p e
      assert(true, "Exception generated and catched")
    end
  end

  #Testst that the coordinate system is valid (can't start from a negative coordinate)
  def test_read_invalid_reference_start_coordinate
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    begin
      als = sam.fetch("chr", -1, 500)
      p als 
      sam.close
      assert(false, "Seems it ran the query")
    rescue Bio::DB::SAMException => e
      p e
      assert(true, "Exception generated and catched")
    end
  end

 #Tests that the query fails when the chromosome doesnt exist and the coordiate is invalid. 
  def test_read_invalid_reference_end_coordinate
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    begin
      als = sam.fetch("chr", 0, 50000)
      p als 
      sam.close
      assert(false, "Seems it ran the query")
    rescue  Bio::DB::SAMException => e
      p e
      assert(true, "Exception generated and catched")
    end
  end
  
  #Tests that an exception is thrown whenver the coordinates are not from left to right. 
  def test_read_invalid_reference_swaped_coordinates
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    begin
      als = sam.fetch("chr", 500, 0)
      p als 
      sam.close
      assert(false, "Seems it ran the query")
    rescue  Bio::DB::SAMException => e
      p e
      assert(true, "Exception generated and catched")
    end
  end

  #Test the load of the fasta index. 
  def test_fasta_load_index
    sam = Bio::DB::Sam.new({:fasta=>@testReference})
    sam.load_reference
    seq = sam.fetch_reference("chr_1", 0, 500)
    p seq 
    sam.close
    assert(true, "The reference was loaded")
  end
  
  #Test that an exeption is thrown when the chromosome is invalid 
  def test_fasta_load_invalid_index
    sam = Bio::DB::Sam.new({:fasta=>@testReference})
    sam.load_reference
    begin
      seq = sam.fetch_reference("chr1", 0, 500)
      p "Error seq:"+ seq 
      sam.close
      assert(false, "The reference was loaded")
    rescue Bio::DB::SAMException => e
      p e
      assert(true,  "The references was not loaded")
    end
  end
  
  #Tests that a feature is loaded, by mapping to a new object outside the library
  def test_load_feature
    fs = Feature.find_by_bam("chr_1", 0, 500,@testBAMFile)
    p fs
    assert(true, "Loaded as features")
  end
  
  #Tests the average coverage function
  def test_avg_coverage
    sam = Bio::DB::Sam.new({:fasta=>@testReference, :bam=>@testBAMFile })
    sam.open
    cov = sam.average_coverage("chr_1", 60, 30)
    p "Coverage: " + cov.to_s
    sam.close
    assert(true, "Average coverage ran")
    assert(3 == cov, "The coverage is 3")
  end
  
  #Test that the average coverage is run with a prefetched index. 
  def test_avg_coverage_cached_name
    sam = Bio::DB::Sam.new({:fasta=>@testReference, :bam=>@testBAMFile })
    sam.open
    sam.fetch_chr_index("chr_1")
  
    p "The index for chr_1: " + ind.to_s
    cov = sam.average_coverage("chr_1", 60, 30)
    p "Coverage: " + cov.to_s
    sam.close
    assert(true, "Average coverage ran")
    assert(3 == cov, "The coverage is 3")
  end

  #Thests the coverage for each possition in the coverage
  def test_chromosome_coverage
    sam = Bio::DB::Sam.new({:fasta=>@testReference, :bam=>@testBAMFile })
    sam.open
    covs = sam.chromosome_coverage("chr_1", 0, 60)
    p "Coverage: "
    p covs
    puts "POS\tCOV"
    covs.each_with_index{ |cov, i| puts "#{i}\t#{cov}" }
    sam.close
    assert(true, "Chromosome coverage ran")
    #assert(3 == cov, "The coverage is 3")
  end
  
  
  
  #Test the sortinf BAM function. 
  def test_sort
     sam = Bio::DB::Sam.new({:fasta=>@testReference, :bam=>@testBAMFile })
     sam.sort(@testBAMFile + "_sorted", 0)
     sam.sort(@testBAMFile + "_sorted_by_name", 1)
     assert(true, "Sorted bam file test ran.")
  end
  #Test that merging works (we still have an excption. I think it works with the latest samtools library, but i have to recompile it with debugging options to be sure)
  def test_merge
     #(files, merged_file, headers, add_RG, by_qname)
     files = []
     files[0] = @testBAMFile
     files[1] = @testBAMFile
     #Bio::DB::Sam.merge(files, "merged.bam", files[0], 0, 0)
     #assert(true, "Sorted bam file test ran.")
  end

end

#Test simplified feature to test that you can transform fromt he default alignment objects to a different mapping
class Feature 
attr_reader :start, :end, :strand, :sequence, :quality

def initialize(a={})
  p a
  @start = a[:start]
  @end = a[:enf]
  @strand = a[:strand]
  @sequence = a[:sequence]
  @quality = a[:quality]
end

def self.find_by_bam(reference,start,stop,bam_file_path)
  
  sam = Bio::DB::Sam.new({:bam=>bam_file_path})
  features = []
  sam.open
  
  fetchAlignment = Proc.new do |a|
    a.query_strand ? strand = '+'  : strand = '-'
    features << Feature.new({:start=>a.pos,:end=>a.calend,:strand=>strand,:sequence=>a.seq,:quality=>a.qual})
  end
  sam.fetch_with_function(reference, start, stop, fetchAlignment)
  sam.close
  features
end
end
