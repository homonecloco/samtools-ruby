# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{samtools-ruby}
  s.version = "0.1"

  s.required_rubygems_version = Gem::Requirement.new(">= 1.2") if s.respond_to? :required_rubygems_version=
  s.authors = ["Ricardo H. Ramirez-Gonzalez"]
  s.date = %q{2010-08-09}
  s.description = %q{Binder of samtools for ruby, on the top of FFI. }
  s.email = %q{}
  s.extra_rdoc_files = ["LICENSE.txt", "README", "lib/bio/db/sam.rb", "lib/bio/db/sam/bam.rb", "lib/bio/db/sam/faidx.rb", "lib/bio/db/sam/sam.rb", "libbam.dylib"]
  s.files = ["FFIGUESS/bam.rb", "FFIGUESS/faidx.rb", "FFIGUESS/sam.rb", "LICENSE.txt", "README", "Rakefile", "lib/bio/db/sam.rb", "lib/bio/db/sam/bam.rb", "lib/bio/db/sam/faidx.rb", "lib/bio/db/sam/sam.rb", "libbam.dylib", "test/basictest.rb", "test/feature.rb", "test/samples/small/sorted.bam", "test/samples/small/test", "test/samples/small/test.bam", "test/samples/small/test.bam.gz", "test/samples/small/test.fa", "test/samples/small/test.fai", "test/samples/small/test.sai", "test/samples/small/test.tam", "test/samples/small/test_chr.fasta", "test/samples/small/test_chr.fasta.amb", "test/samples/small/test_chr.fasta.ann", "test/samples/small/test_chr.fasta.bwt", "test/samples/small/test_chr.fasta.pac", "test/samples/small/test_chr.fasta.rbwt", "test/samples/small/test_chr.fasta.rpac", "test/samples/small/test_chr.fasta.rsa", "test/samples/small/test_chr.fasta.sa", "test/samples/small/testu.bam", "test/samples/small/testu.bam.bai", "Manifest", "samtools-ruby.gemspec"]
  s.homepage = %q{http://github.com/homonecloco/samtools-ruby}
  s.rdoc_options = ["--line-numbers", "--inline-source", "--title", "Samtools-ruby", "--main", "README"]
  s.require_paths = ["lib"]
  s.rubyforge_project = %q{samtools-ruby}
  s.rubygems_version = %q{1.3.6}
  s.summary = %q{Binder of samtools for ruby, on the top of FFI.}

  if s.respond_to? :specification_version then
    current_version = Gem::Specification::CURRENT_SPECIFICATION_VERSION
    s.specification_version = 3

    if Gem::Version.new(Gem::RubyGemsVersion) >= Gem::Version.new('1.2.0') then
    else
    end
  else
  end
end
