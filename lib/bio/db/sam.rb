require 'bio/db/sam/sam'
module Bio
  module DB
    class Sam
      attr_reader :sam_file

      def initialize(optsa={})
        opts =  { :fasta => nil,  :bam => nil,:tam => nil, :compressed => true, :write => false }.merge!(optsa)



        @fasta_path = opts[:fasta]
        @compressed = opts[:compressed]
        @write = opts[:write]
        bam = opts[:bam]
        tam = opts[:tam]

        if bam == nil && tam == nil then
          raise SAMException.new(), "No alignment file"
        elsif bam != nil && tam != nil then
          raise SAMException.new(), "Alignment has to be in either text or binary format, not both"
        elsif bam != nil then
          @binary = true
          @sam = bam     
        else
          @sam = tam     
          @binary = false
        end
        @sam_file = nil

        ObjectSpace.define_finalizer(self,  self.class.method(:finalize).to_proc)
      end

      def open()
        if(@write)
          raise SAMException.new(), "Writing not supported"
        end

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

      def to_s()
        (@binary ? "Binary" : "Text") + " file: " + @sam + " with fasta: " + @fasta_path
      end

      def close()

        Bio::DB::SAM::Tools.bam_index_destroy(@sam_index) unless @sam_index.nil? || @sam_index.null?
        Bio::DB::SAM::Tools.samclose(@sam_file) unless @sam_file.nil? 
        @sam_file = nil
      end

      def Sam.finalize(id)
        id.close()
        puts "Finalizing #{id}  at #{Time.new}"       
      end

      def load_index()
        raise SAMException.new(), "Indexes are only supported by BAM files, please use samtools to convert your SAM file" unless @binary
        @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
        if @sam_index.null? then
          p "Generatind index for: " + @sam
          Bio::DB::SAM::Tools.bam_index_build(@sam)
          @sam_index = Bio::DB::SAM::Tools.bam_index_load(@sam)
          raise SAMException.new(), "Unable to generate bam index for: " + @sam if @sam_index.null?
        end
      end

      
      
      def fetch(chromosome, qstart, qend)
        load_index if @sam_index.nil? || @sam_index.null?
        chr = FFI::MemoryPointer.new :int
        beg = FFI::MemoryPointer.new :int
        last = FFI::MemoryPointer.new :int
        query =   chromosome+":"+qstart.to_s+"-"+qend.to_s 
        qpointer = FFI::MemoryPointer.from_string(query)
        header = @sam_file[:header]
        Bio::DB::SAM::Tools.bam_parse_region(header,qpointer, chr, beg, last) 
        als = Array.new
        fetchAlignment = Proc.new do |bam_alignment, data|
          alignment =  Alignment.new
          alignment.set(bam_alignment, header)
          als.push(alignment)   
          als 
	  0  
        end
        Bio::DB::SAM::Tools.bam_fetch(@sam_file[:x][:bam], @sam_index,chr.read_int,beg.read_int, last.read_int, nil, fetchAlignment)
       # p als
        #puts "Iterated"
        als
      end    

    end
    class Tag
      attr_accessor :tag, :type, :value
      def set(str)
        v = str.split(":")
        @tag   = v[0]
        @type  = v[1]
        @value = v[2]
      end
    end
    
    class Alignment
       attr_accessor :qname, :flag, :rname,:pos,:mapq,:cigar, :mrnm, :mpos, :isize, :seq, :qual, :tags, :calend, :qlen
       
       def set(bam_alignment, header)
         al = Bio::DB::SAM::Tools::Bam1T.new(bam_alignment) 
         self.sam =  Bio::DB::SAM::Tools.bam_format1(header,al)
         core = al[:core]
          cigar = al[:data][core[:l_qname]]#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
      
          
         @calend = Bio::DB::SAM::Tools.bam_calend(core,cigar)
         @qlen = Bio::DB::SAM::Tools.bam_cigar2qlen(core,cigar)
         
       end
       
    #   def sam
     #      @sam
      # end
       
       def sam=(sam)
           p sam
           s = sam.split("\t")
           self.qname = s[0]
           self.flag  = s[1]
           self.rname = s[2]
           self.pos   = s[3]
           self.mapq  = s[4]
           self.cigar = s[5]
           self.mrnm  = s[6]
           self.mpos  = s[7]
           self.isize = s[8]
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

    class SAMException < RuntimeError
      #we can add further variables to give information of the excpetion
      def initialize()

      end
    end
  end
end


