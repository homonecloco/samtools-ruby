require 'rubygems'
require'ffi'
require 'bio/db/sam/bam'
# This module is a direct mapping of sam.h from the samtools library
module Bio
  module DB
    module SAM
      module Tools
        extend FFI::Library
        ffi_lib'libbam'
        class BamHeaderT < FFI::Struct
          layout(
          :n_targets, :int32_t,
          :target_name, :pointer,
          :target_len, :pointer,
          :dict, :pointer,
          :hash, :pointer,
          :rg2lib, :pointer,
          :l_text, :int,
          :text, :pointer
          )
          def text=(str)
            @text = FFI::MemoryPointer.from_string(str)
            self[:text] = @text
          end
          def text
            @text.get_string(0)
          end

        end
        
        class SamfileTX < FFI::Union
          layout(
          :tamr, :pointer, #Text file, read.
          :bam,  :pointer, #bamFile,
          :tamw, :pointer #Text file, write. 
          )
        end
        class SamfileT < FFI::Struct
          layout(
          :type, :int,          
          :x, SamfileTX,
          :header, :pointer
          )
        end
        
        

        attach_function :samclose, [ :pointer ], :void
        attach_function :samread, [ :pointer, :pointer ], :int
        attach_function :samopen, [ :string, :string, :pointer ], :pointer
        attach_function :samwrite, [ :pointer, :pointer ], :int
        attach_function :sampileup, [ :pointer, :int, :bam_pileup_f, :pointer ], :int
        attach_function :samfaipath, [ :string ], :string
      end
    end
  end
end

