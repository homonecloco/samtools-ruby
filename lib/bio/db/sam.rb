require 'bio/db/sam/sam'
require 'bio/db/sam/libc'


#This module is a high level binder for the libbam library from samtools. 
#The code is written to hide most of the low level complexity and give a more ruby-like API. 

module Bio
  module DB
    class Sam
      attr_reader :sam_file, :chr_index
      attr_accessor :sort_mem_size

      #To make a new sam object. Initialize expects a hash optsa with the following elemets:
      # fasta:: The fasta file with the reference. (nil)
      # bam:: path to a binary SAM file (nil)
      # tam:: path to a text SAM file (nil) 
      # compressed:: If the binary file is compressed (true)
      # write:: If the file is to be writen (false). Not supported yet. 
      # *NOTE:* you can't use binary and text formats simultaneusly. To make queries, the file has to be a sorted binary. 
      # This function doesn't actually open the file, it just prepares the object to be opened in a later stage. 
      # 
      def initialize(optsa={})
        opts =  { :fasta => nil,  :bam => nil,:tam => nil, :compressed => true, :write => false }.merge!(optsa)

        @sort_mem_size = 500000000
        @chr_index = {}
        @fasta_path = opts[:fasta]
        @compressed = opts[:compressed]
        @write = opts[:write]
        bam = opts[:bam]
        tam = opts[:tam]

        if bam == nil && tam == nil && @fasta_path == nil then
          raise SAMException.new(), "No alignment or reference file"
        elsif bam != nil && tam != nil then
          raise SAMException.new(), "Alignment has to be in either text or binary format, not both"
        elsif bam != nil then
          @binary = true
          @sam = bam     
        elsif tam != nil then
          @sam = tam     
          @binary = false

        end
        @fasta_file = nil
        @sam_file   = nil

        ObjectSpace.define_finalizer(self,  self.class.method(:finalize).to_proc)
      end
      
      #Function that actually opens the sam file
      #Throws a SAMException if the file can't be open. 
      def open()

        raise SAMException.new(), "Writing not supported yet" if @write
        raise SAMException.new(), "No SAM file specified" unless @sam 

        opts = @write ? "w" : "r"
        if @binary then  
          opts += "b" 
          if @write then
            unless @compressed then 
              opts += "u"
            end
          end
        end
        valid = ["r", "w", "wh", "rb", "wb" , "wbu"]
        unless valid.include?(opts) then
          raise SAMException.new(), "Invalid options for samopen: " + opts 
        end

        samFile = Bio::DB::SAM::Tools.samopen(@sam, opts, nil)
        if samFile.null? then
          @sam_file = nil
          raise SAMException.new(), "File not opened:  " + @sam
        end
        @sam_file = Bio::DB::SAM::Tools::SamfileT.new(samFile)

      end

      #Prints a description of the sam file in a text format containg if it is binary or text, the path
      #and the fasta file of the reference
      def to_s()
        (@binary ? "Binary" : "Text") + " file: " + @sam + " with fasta: " + @fasta_path
      end

      #Closes the sam file and destroys the C pointers using the functions provided by libbam
      def close()
        Bio::DB::SAM::Tools.fai_destroy(@fasta_index) unless @fasta_index.nil? || @fasta_index.null?
        Bio::DB::SAM::Tools.bam_index_destroy(@sam_index) unless @sam_index.nil? || @sam_index.null?
        Bio::DB::SAM::Tools.samclose(@sam_file) unless @sam_file.nil? 
        @sam_file = nil
        @fasta_index = nil
        @chr_index = nil
      end

      # Destructor method that closes the file before letting the object be garbage collected. 
      def Sam.finalize(id)
        id.close()      
      end
      
      #Loads the bam index to be used for fetching. If the index doesn't exists the index is built provided that
      #the user has writing access to the folder where the BAM file is located. If the creation of the file fails
      #a SAMException is thrown. 
      #If the index doesn't exist, loading it will take more time. It is suggested to generate the index separatedly
      #if the bam file sits on a server where the executing user may not have writing permissions in the server.
      def load_index()
        raise SAMException.new(), "Indexes are only supported by BAM files, please use samtools to convert your SAM file" unless @binary
        @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
        if @sam_index.null? then
          p "Generating index for: " + @sam
          Bio::DB::SAM::Tools.bam_index_build(@sam)
          @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
          raise SAMException.new(), "Unable to generate bam index for: " + @sam if @sam_index.nil? || @sam_index.null?
        end
      end

      #Loads the reference file to be able to query regions of it. This requires the fai index to exist in the same
      #folder than the reference. If it doesn't exisits, this functions attempts to generate it. If user doesn't
      #have writing permissions on the folder, or the creation of the fai fails for any reason, a SAMException is thrown.
      def load_reference()
        raise SAMException.new(), "No path for the refernce fasta file. " if @fasta_path.nil?

        @fasta_index = Bio::DB::SAM::Tools.fai_load(@fasta_path)

        if @fasta_index.null? then
          p "Generating index for: " + @fasta_path
          Bio::DB::SAM::Tools.fai_build(@fasta_path)
          @fasta_index =  Bio::DB::SAM::Tools.fai_load(@fasta_path)
          raise SAMException.new(), "Unable to generate fasta index for: " + @fasta_path if @fasta_index.nil? ||  @fasta_index.null?
        end

      end

      #Returns the average coverage of a region in a bam file. 
      def average_coverage(chromosome, qstart, len)

        coverages = chromosome_coverage(chromosome, qstart, len)
        total = 0
        len.times{ |i| total= total + coverages[i] }
        avg_cov = total.to_f / len
        
        avg_cov
      end

      #Returns an array with the coverage at each possition in the queried region
      #This is a simple average coverage just calculated with the first and last
      #possition of the alignment, ignoring the gaps. 
      def chromosome_coverage(chromosome, qstart, len)
        coverages = Array.new(len, 0)
        chr_cov_proc = Proc.new do |alignment|
          last = alignment.calend - qstart
          first = alignment.pos - qstart
          if last < first
            tmp = last
            last = first
            first = last
          end
          first.upto(last-1) { |i|
            coverages[i-1] = 1 + coverages[i-1]  if i-1 < len && i > 0
          }
        end
        fetch_with_function(chromosome, qstart, qstart+len,  chr_cov_proc)
        coverages
      end

      #Returns the sequence for a given region. 
      def fetch_reference(chromosome, qstart,qend)
        load_reference if @fasta_index.nil? || @fasta_index.null? 
        query = query_string(chromosome, qstart,qend)
        len = FFI::MemoryPointer.new :int
        reference = Bio::DB::SAM::Tools.fai_fetch(@fasta_index, query, len)
        raise SAMException.new(), "Unable to get sequence for reference: "+query if reference.nil?
        reference
      end

      #Generates a query sting to be used by the region parser in samtools. 
      #In principle, you shouldn't need to use this function. 
      def query_string(chromosome, qstart,qend)
        query = chromosome + ":" + qstart.to_s + "-" + qend.to_s 
        query
      end

      #Returns an array of Alignments on a given region. 
      def fetch(chromosome, qstart, qend)
        als = Array.new
        fetchAlignment = Proc.new do |alignment|
          als.push(alignment.clone)   
          0  
        end
        fetch_with_function(chromosome, qstart, qend, fetchAlignment)
        als
      end  


      #Gets the internal index for a chromosome in the reference. 
      #This function calls the C API to get a pointer to the index. 
      #Throws a SAMException if the chromosme is not in the index.  
      #This shouldn't be done outside this Module.
      def fetch_chr_index(chromosome)
        chr = FFI::MemoryPointer.new :int
        beg = FFI::MemoryPointer.new :int
        last = FFI::MemoryPointer.new :int
        
        header = @sam_file[:header]
        
        query = query_string(chromosome, 1,2)
        qpointer = FFI::MemoryPointer.from_string(query)
        Bio::DB::SAM::Tools.bam_parse_region(header,qpointer, chr, beg, last) 
        raise SAMException.new(), "invalid chromosome: " + chromosome  if(chr.read_int < 0)
        @chr_index[chromosome] = chr.read_int
        @chr_index[chromosome]
      end

      #Gets the internal index for a chromosome. This function has a caching mechanism
      #to avoid quering each time for the index, so the C API is called only the first
      #time a chromosome index is required.
      #This shouldn't be done outside this Module.
      def get_chr_index(chromosome)
         load_index if @sam_index.nil? || @sam_index.null?
         fetch_chr_index(chromosome) unless @chr_index[chromosome]
         @chr_index[chromosome]
      end

      #Executes a function on each Alignment inside the queried region of the chromosome. The chromosome
      #can be either the textual name or a FixNum with the internal index. However, you need to get the
      #internal index with the provided API, otherwise the pointer is outside the scope of the C library. 
      #Returns the count of alignments in the region. 
      #WARNING: Accepts an index already parsed by the library. It fails when you use your own FixNum (FFI-bug?)
      def fetch_with_function(chromosome, qstart, qend, function)
        tmp = "chromosome.class: " + chromosome.class.to_s
        p tmp
        chr_i = chromosome if chromosome.instance_of?(Fixnum)
        chr_i = get_chr_index(chromosome) unless chromosome.instance_of?(Fixnum)
        header = @sam_file[:header]
        
        count = 0;

        fetchAlignment = Proc.new do |bam_alignment, data|
          alignment =  Alignment.new
          alignment.set(bam_alignment, header)
          function.call(alignment)
          count = count + 1
          0  
        end
        Bio::DB::SAM::Tools.bam_fetch(@sam_file[:x][:bam], @sam_index,chr_i,qstart, qend, nil, fetchAlignment)
        count
      end
      
      #Sorts the sam file and stores it in a new file.
      #prefix:: the prefix of the sorted file
      #by_qname:: if true, the file is sorted by query name. 
      def sort(prefix, by_qname)
         #void bam_sort_core_ext(int void bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int add_RG), const char *fn, const char *prefix, size_t max_mem, int is_stdout)
        Bio::DB::SAM::Tools.bam_sort_core_ext(by_qname, @sam, prefix, @sort_mem_size, 0)
      end
      
      #Merges n BAM files. This doesn't require to create a SAM object
      #files:: An array with the paths to the files.
      #merged_file:: The path to the merged file
      #headers:: The BAM file containing the header
      #add_RG:: If true, the RG tag is added (infered from the filenames)
      #by_qname:: If true, the bamfiles should by ordered by query name, if false, by coordinates. 
      def self.merge(files, merged_file, headers, add_RG, by_qname)
        strptrs = []
        files.each do |file|
          strptrs << FFI::MemoryPointer.from_string(file)
        end
        strptrs << nil

        # Now load all the pointers into a native memory block
        argv = FFI::MemoryPointer.new(:pointer, strptrs.length)
        strptrs.each_with_index do |p, i|
           argv[i].put_pointer(0,  p)
        end
        #void bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int add_RG)
        Bio::DB::SAM::Tools.bam_merge_core(by_qname, merged_file, headers, strptrs.length, argv, add_RG)
      end
    end

    #Contains the tags for the alignments. See the SAM format documentation. 
    class Tag
      attr_accessor :tag, :type, :value
      def set(str)
        v = str.split(":")
        @tag   = v[0]
        @type  = v[1]
        @value = v[2]
      end
    end

    #Contains a single alignment. It wrapps all the infomration defined by a single Alignment in the SAM format. 
    #Attrobhtes frp, the flag field (see chapter 2.2.2 of the sam file documentation)
    #query_strand and mate_strand are true if they are forward. It is the opposite to the definition in the BAM format for clarity.
    #primary is the negation of is_negative from the BAM format
    #Besides the stored information, it has the calend and qlen from the API. 
    class Alignment

      #Constructor that sets the destructor
      def initialize
        ObjectSpace.define_finalizer(self,
        self.class.method(:finalize).to_proc)
      end
      
      #This free the temporary memory used by the alignment
      def Alignment.finalize(object_id)
        LibC.free object_id.al
        LibC.free object_id.sam
        LibC.free object_id.calend
        LibC.free object_id.qlen
        LibC.free object_id.samstr
      end

      #Attributes from the format
      attr_accessor :qname, :flag, :rname,:pos,:mapq,:cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags, :al, :samstr
      #Attributes pulled with the C library
      attr_accessor  :calend, :qlen
      
      attr_accessor :is_paired, :is_mapped, :query_unmapped, :mate_unmapped, :query_strand, :mate_strand, :first_in_pair,:second_in_pair, :primary, :failed_quality, :is_duplicate

      #Sets the alignment from a pointer to BAM alignment from the C library. 
      def set(bam_alignment, header)
        #Create the FFI object
        @al = Bio::DB::SAM::Tools::Bam1T.new(bam_alignment) 

        #set the raw data
        tmp_str =  Bio::DB::SAM::Tools.bam_format1(header,al)
      
        self.sam = String.new(tmp_str)
        
        #Set values calculated by libbam
        core = al[:core]
        cigar = al[:data][core[:l_qname]]#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname)) 
        @calend = Bio::DB::SAM::Tools.bam_calend(core,cigar)
        @qlen = Bio::DB::SAM::Tools.bam_cigar2qlen(core,cigar)

        #process the flags
        @is_paired             = @flag & 0x0001 > 0
        @is_mapped             = @flag & 0x0002 > 0
        @query_unmapped        = @flag & 0x0004 > 0
        @mate_unmapped         = @flag & 0x0008 > 0
        @query_strand          = !(@flag & 0x0010 > 0)
        @mate_strand           = !(@flag & 0x0020 > 0)
        @first_in_pair         = @flag & 0x0040 > 0
        @second_in_pair        = @flag & 0x0080 > 0
        @primary               = !(@flag & 0x0100 > 0)
        @failed_quality        = @flag & 0x0200 > 0
        @is_duplicate          = @flag & 0x0400 > 0

      end

      #Sets the sam values from a SAM line. 
      def sam=(sam)
        #p sam
        s = sam.split("\t")
        self.qname = s[0]
        self.flag  = s[1].to_i
        self.rname = s[2]
        self.pos   = s[3].to_i
        self.mapq  = s[4].to_i
        self.cigar = s[5]
        self.mrnm  = s[6]
        self.mpos  = s[7].to_i
        self.isize = s[8].to_i
        self.seq   = s[9]
        self.qual =  s[10]
        self.tags = {} 
        11.upto(s.size-1) {|n| 
          t = Tag.new 
          t.set(s[n])
          tags[t.tag] = t
        }


        #<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
        #[<TAG>:<VTYPE>:<VALUE> [...]] 

      end

    end

    #Exception thrown whenver something fails inside the API. 
    class SAMException < RuntimeError
      #we can add further variables to give information of the excpetion
      def initialize()

      end
    end
  end
end


