require 'bio'

class Bio::Velvet::Graph::NodedRead
  def adjusted_position(parent_node)
    if @direction == true
      return @offset_from_start_of_node
    elsif @direction == false
      return parent_node.length - @offset_from_start_of_node
    else
      raise "programming error"
    end
  end
end

module Bio
  module AssemblyGraphAlgorithms
    class ContigPrinter
      include Bio::FinishM::Logging

      class AnchoredConnection
        # The identifiers of the probe reads in the velvet assembly graph
        attr_accessor :start_probe_noded_read, :end_probe_noded_read

        # number of nucleotides between the start of the start probe read and the start of the end of the contig
        attr_accessor :start_probe_contig_offset

        # number of nucleotides until the end of the end probe read in the start of the second contig
        attr_accessor :end_probe_contig_offset

        # Length of the start and end probe sequences
        attr_accessor :start_probe_read_length
        attr_accessor :end_probe_read_length

        # Enumerable of Enumerables of OrientedNode objects, each list of OrientedNode objects
        # corresponds to a path that forms the connection
        attr_accessor :paths
      end

      # Much like one_connection_between_two_contigs except can handle multiple connections
      # (but cannot handle 0 connections)
      def ready_two_contigs_and_connections(graph, contig1, anchored_connection, contig2, sequences)
        to_return = ''
        variants = []

        log.debug "Working with anchored_connection: #{anchored_connection.inspect}" if log.debug?

        # Stage1 - contig1 before the path begins
        to_return = nil
        if anchored_connection.start_probe_contig_offset == 0
          # 0 is a special case because negative 0 doesn't make sense
          to_return = contig1
        else
          to_return = contig1[0...-(anchored_connection.start_probe_contig_offset)]
        end
        log.debug "After first chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?

        # Stage2 - path sequence, beginning and ending with
        # beginning and ending probes
        begin
          reference_path_index = predict_reference_path_index(anchored_connection.paths)
          path = anchored_connection.paths[reference_path_index]
          path_sequence = path.sequence
          log.debug "Reference path has a sequence length #{path_sequence.length}" if log.debug?

          # Find start index
          begin_onode = path[0]
          begin_noded_read = anchored_connection.start_probe_noded_read
          raise if begin_noded_read.nil?
          extra_bit_on_start = ''
          if begin_noded_read.start_coord != 0
            log.warn "Unexpectedly the start of the start probe not did not form part of the path, which is a little suspicious"
            extra_bit_on_start = sequences[begin_noded_read.read_id][0...begin_noded_read.start_coord]
          end
          offset_of_begin_probe_on_path = nil
          # xor read direction on node, and node direction on path
          if (begin_noded_read.direction == true) ^ begin_onode.starts_at_start?
            offset_of_begin_probe_on_path = begin_onode.node.corresponding_contig_length - begin_noded_read.offset_from_start_of_node
            # extra bit on read needs to be reverse complemented
            extra_bit_on_start = Bio::Sequence::NA.new(extra_bit_on_start).reverse_complement.to_s.upcase unless extra_bit_on_start == ''
          else
            offset_of_begin_probe_on_path = begin_noded_read.offset_from_start_of_node
          end

          # Work out the variants if there is any
          variants = paths_to_variants(
            path,
            anchored_connection.paths - [path] #Ain't Ruby grand..
            )
          # Correct variants' positions to be relative to the full contig,
          # not just the path sequence
          variants.each do |variant|
            variant.position = variant.position - offset_of_begin_probe_on_path + to_return.length + 1
          end

          # Find end index
          end_onode = path[-1]
          end_noded_read = anchored_connection.end_probe_noded_read
          raise if end_noded_read.nil?
          extra_bit_on_end = ''
          if end_noded_read.start_coord != 0
            log.warn "Unexpectedly the end of the end probe not did not form part of the path, which is a little suspicious"
            extra_bit_on_end = sequences[end_noded_read.read_id][0...end_noded_read.start_coord]
          end
          offset_of_end_node_on_path = path[0...-1].reduce(0){|sum, onode| sum += onode.node.length_alone}
          if (end_noded_read.direction == false) ^ end_onode.starts_at_start?
            offset_of_end_node_on_path += end_noded_read.offset_from_start_of_node
            extra_bit_on_end = Bio::Sequence::NA.new(extra_bit_on_end).reverse_complement.to_s.upcase unless extra_bit_on_end == ''
          else
            offset_of_end_node_on_path += end_onode.node.corresponding_contig_length - end_noded_read.offset_from_start_of_node
          end

          log.debug "Found start index #{offset_of_begin_probe_on_path} and end index #{offset_of_end_node_on_path}" if log.debug?
          to_return += extra_bit_on_start+
            path_sequence[offset_of_begin_probe_on_path...offset_of_end_node_on_path]+
            extra_bit_on_end
          log.debug "After path chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?
        end #end stage 2

        # Stage 3
        to_return += contig2[anchored_connection.end_probe_contig_offset..-1]
        log.debug "After last chunk of sequence added, sequence is #{to_return.length}bp long" if log.debug?

        return to_return, variants
      end

      # Given two contigs, return a String representing the new contig. Assumes
      # that there is only 1 path between the two contigs.
      #
      #          ---------->         <--------           start and end probes (ends of probe sequences may not form part of final path). Directions not variable.
      #  --------------------->NNNN------------------->  original sequence to be gapfilled (contig1, NNNN, contig2). Directions not variable
      #      -----------                 ------->        path across the gap. Direction not variable
      #                 \               /
      #                  --------------
      #      ---------->|<-----|----->|--------->        nodes that make up the path (directions and boundaries variable)
      #    stage1|           stage2           |stage3    stages of sequence construction in this method
      def one_connection_between_two_contigs(graph, contig1, anchored_connection, contig2, sequences)
        raise "programming error: only one path expected here" if anchored_connection.paths.length > 1
        return ready_two_contigs_and_connections(graph, contig1, anchored_connection, contig2, sequences)[0]
      end

      # Return the index of a reference path picked from the given paths.
      # Current method is simply highest coverage path
      def predict_reference_path_index(paths)
        max_i = 0
        max_coverage = paths[0].coverage
        paths[1..-1].each_with_index do |path, i|
          cov = path.coverage
          if cov > max_coverage
            max_i = i+1
            max_coverage = cov
          end
        end
        return max_i
      end

      # Return an Array of Variant objects
      def paths_to_variants(reference_path, non_reference_paths)
        sequences_to_variants(
          reference_path.sequence,
          non_reference_paths.collect{|path| path.sequence}
          )
      end

      def sequences_to_variants(reference_sequence, alternate_sequences)
        return [] if alternate_sequences.empty?

        # Run multiple sequence alignment of each sequence, with the reference sequence first
        log.debug "Running MSA with #{1+alternate_sequences.length} sequences.." if log.debug?
        alignments = clustalo([
          reference_sequence,
          alternate_sequences
          ].flatten)
        log.debug "Finished running MSA" if log.debug?
        if log.debug?
          log.debug "Alignment found was:"
          alignments.each do |align|
            log.debug align
          end
        end

        # Collect the variants at each sequence at each column
        ref_alignment = alignments[0]
        non_ref_alignments = alignments[1..-1]
        variants = [] #Array of empty arrays
        reference_position = 0
        i = 0
        ref_alignment.each_char do |ref_base|
          non_ref_alignments.each_with_index do |alignment, sequence_id|
            nonref = alignment[i]
            if nonref != ref_base
              variant = nil
              if ref_base == '-'
                variant = Variant.new reference_position, nonref, Variant::INSERT
              elsif nonref == '-'
                variant = Variant.new reference_position, 1, Variant::DELETION
              else
                variant = Variant.new reference_position, nonref, Variant::SWAP
              end
              variants[sequence_id] ||= []
              variants[sequence_id].push variant
            end
          end
          reference_position += 1 unless ref_base == '-'
          i += 1
        end
        #        if log.debug?
        #          log.debug "Before variation, found #{variants.length} variants:"
        #          variants.each_with_index do |variant_set, alt_i|
        #            variant_set.each do |variant|
        #              log.debug "From #{alt_i}: #{variant.to_shorthand}"
        #            end
        #          end
        #        end

        # Condense the single column, single species variants into a condensed set
        return condense_variants!(variants)
      end

      def condense_variants!(variant_array_of_arrays)
        all_variants = {}

        variant_array_of_arrays.each_with_index do |variant_array, i|
          last_variant = nil
          current_variants = []
          variant_array.each do |variant|
            # Combine last_variant and this one if
            # their positions are consecutive and their types are the same
            if !last_variant.nil? and last_variant.type == variant.type

              if variant.type == Variant::INSERT and last_variant.position == variant.position
                last_variant.sequence += variant.sequence

              elsif variant.type == Variant::DELETION and last_variant.position == variant.position - last_variant.deletion_length
                last_variant.deletion_length += 1

              elsif variant.type == Variant::SWAP and last_variant.position + last_variant.sequence.length == variant.position
                last_variant.sequence += variant.sequence

              else
                # Start a new variant
                last_variant = variant
                current_variants.push variant
              end
            else
              last_variant = variant
              current_variants.push variant
            end
          end
          if log.debug?
            log.debug "Found #{current_variants.length} variants in sequence #{i}:"
            current_variants.each do |variant|
              log.debug variant.to_shorthand
            end
          end

          # Multiple paths can have the same variant. Don't duplicate
          current_variants.each do |variant|
            key = [
              variant.position,
              variant.sequence,
              variant.deletion_length,
              variant.type
              ]
            all_variants[key] ||= variant
          end
        end

        return all_variants.values
      end

      #       # Given an Enumerable of nucleic acid sequences, align them with MAFFT,
      #       # and return an Array of the same size as the input
      #       def mafft(sequences)
      #         i = 0
      #         stdin = sequences.collect{|s| i+=1; ">#{i}\n#{s}\n"}.join('')
      #         stdout = Bio::Commandeer.run "mafft --retree 1 --quiet --nuc /dev/stdin", {:stdin => stdin, :log => log}
      #         to_return = []
      #         header = true
      #         stdout.each_line do |line|
      #           if !header
      #             to_return.push line.strip
      #           end
      #           header = !header
      #         end
      #         return to_return
      #       end

      def clustalo(sequences)
        i = 0
        stdin = sequences.collect{|s| i+=1; ">#{i}\n#{s}\n"}.join('')
        stdout = Bio::Commandeer.run "clustalo -t DNA -i - --output-order=input-order", {:stdin => stdin, :log => log}
        to_return = []
        header = true
        Bio::FlatFile.foreach(Bio::FastaFormat, StringIO.new(stdout)) do |seq|
          to_return.push seq.seq.to_s
        end
        return to_return
      end




      class Variant
        #Types:
        INSERT = :insert
        DELETION = :deletion
        SWAP = :swap #n bases swapped for another n bases

        attr_accessor :reference_name

        # 0-based position on the contig
        attr_accessor :position

        # sequence (or nil if variant is a deletion)
        attr_accessor :sequence

        # length of deletion (or nil if not a deletion)
        attr_accessor :deletion_length

        # See constants in this class
        attr_accessor :type

        def initialize(position=nil, sequence_or_deletion_length=nil, type=nil)
          @position = position
          @type = type
          if type == DELETION
            @deletion_length = sequence_or_deletion_length
          else
            @sequence = sequence_or_deletion_length
          end
        end

        def base_number
          @position+1
        end

        def to_shorthand
          if type == DELETION
            "#{base_number}D:#{deletion_length}"
          elsif type == SWAP
            "#{base_number}S:#{sequence.upcase}"
          elsif type == INSERT
            "#{base_number}I:#{sequence.upcase}"
          else
            raise
          end
        end

        # The reference sequence has been reverse complemented. Fix this
        # variant so it makes sense again (position aside)
        def reverse!
          if type == SWAP or type == INSERT
            @sequence = Bio::Sequence::NA.new(@sequence).reverse_complement.to_s.upcase
          end
        end

        #CHROM POS     ID        REF    ALT     QUAL FILTER INFO
        def vcf_array(reference_sequence)
          bits = [
            @reference_name,
            @position+1,
            '.',
            ]
          case type
          when SWAP then
            bits.push reference_sequence[@position...(@position+@sequence.length) ]
            bits.push @sequence
          when INSERT then
            bits.push '.'
            bits.push @sequence
          when DELETION then
            bits.push reference_sequence[@position...(@position+@deletion_length) ]
            bits.push '.'
          else
            raise
          end

            bits.push '20'
            bits.push 'PASS'
            bits.push 'finishm'
            return bits
        end

            def vcf(reference_sequence)
              vcf_array(reference_sequence).join("\t")
            end
      end

      class PrintableConnection
        attr_accessor :reference_path, :variants
      end
    end
  end
end
