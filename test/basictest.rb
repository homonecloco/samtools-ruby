$: << File.expand_path(File.dirname(__FILE__) + '/../lib')
require "test/unit"
require "bio/db/sam"
require "bio/db/sam/sam"

class TestBioDbSam < Test::Unit::TestCase


  def setup
    @testTAMFile                 = "test/samples/small/test.tam"
    @testBAMFile                 = "test/samples/small/testu.bam"
  end

  def default_test
    puts $LOAD_PATH
    assert(true, "Unit test test")
  end

  def test_openSAMFile
    bamfile                      = Bio::DB::SAM::Tools.samopen(@testTAMFile,"r",nil)
    Bio::DB::SAM::Tools.samclose(bamfile)
    assert(true, "file open and closed")
  end

  def test_new_class_empty
    begin
      bam                        = Bio::DB::Sam.new({})
      assert(false, "Should fail while opening without parameters")
    rescue Bio::DB::SAMException => e
      puts e.message
      assert(true, e.message)
    end
  end

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

  def test_class_text_read_no_faidx
    sam                          = Bio::DB::Sam.new({:tam=>@testTAMFile})
    sam.open
    sam.close
    assert(true, "file open and closed with the class")
  end

  def test_class_text_read_no_close

    fam                          = Bio::DB::Sam.new({:tam=>@testTAMFile})
    fam.open
    fam                          = nil   
    ObjectSpace.garbage_collect

    assert(true, "file openend but not closed")
  end

  def test_class_binary_read_no_close

    Bio::DB::Sam.new({:bam=>@testBAMFile}).open
    ObjectSpace.garbage_collect
    assert(true, "BINARY file openend but not closed")  
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

  def test_bam_load_index
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    index = sam.load_index
    sam.close
    assert(true, "BAM index loaded")
    #  attach_function :bam_index_build, [ :string ], :int
    #  attach_function :bam_index_load, [ :string ], :pointer
    #  attach_function :bam_index_destroy, [ :pointer ], :void
  end

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

  def test_read_segment
    sam       = Bio::DB::Sam.new({:bam=>@testBAMFile})
    sam.open
    als = sam. fetch("chr_1", 0, 500)
    sam.close
    assert(true, "Seems it ran the query")
    #node_7263       238     60 has 550+, query from 0 to 500, something shall come.... 
  end

end
