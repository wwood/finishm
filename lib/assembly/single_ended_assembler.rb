require 'ds'
require 'set'
require 'ruby-progressbar'

class Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  DEFAULT_MAX_TIP_LENGTH = 200
  DEFAULT_MIN_CONTIG_SIZE = 500
  DEFAULT_MIN_CONFIRMING_RECOHERENCE_READS = 2

  attr_accessor :graph

  ASSEMBLY_OPTIONS = [
    :max_tip_length,
    :recoherence_kmer,
    :min_confirming_recoherence_kmer_reads,
    :sequences,
    :leash_length,
    :min_contig_size,
    :max_coverage_at_fork,
    ]
  attr_accessor :assembly_options

  # Create a new assembler given a velvet graph and velvet Sequences object
  #
  # Assembly options:
  # :max_tip_length: if a path is shorter than this in bp, then it will be clipped from the path. Default 100
  # :recoherence_kmer: attempt to separate paths by going back to the reads with this larger kmer (requires :seqeunces)
  # :sequences: the sequences of the actual reads, probably a Bio::Velvet::Underground::BinarySequenceStore object
  # :leash_length: don't continue assembly from nodes farther than this distance (in bp) away
  # :min_coverage_of_start_nodes: only start exploring from nodes with this much coverage
  # :min_contig_size: don't bother returning contigs shorter than this (default 500bp)
  # :progressbar_io: given an IO object e.g. $stdout, write progress information
  def initialize(graph, assembly_options={})
    @graph = graph
    @assembly_options = assembly_options
    @assembly_options[:max_tip_length] ||= DEFAULT_MAX_TIP_LENGTH
    @assembly_options[:min_contig_size] ||= DEFAULT_MIN_CONTIG_SIZE
    @assembly_options[:min_confirming_recoherence_kmer_reads] ||= DEFAULT_MIN_CONFIRMING_RECOHERENCE_READS
  end

  # Assemble everything in the graph into OrientedNodeTrail objects.
  # Yields an OrientedNodeTrail if a block is
  # given, otherwise returns an array of found paths. Options for
  # assembly are specified in assembly_options
  def assemble
    paths = []

    # Gather a list of nodes to try starting from
    starting_nodes = gather_starting_nodes
    log.info "Found #{starting_nodes.length} nodes to attempt assembly from"

    seen_nodes = Set.new
    progress = setup_progressbar starting_nodes.length

    # For each starting node, start the assembly process
    dummy_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    starting_nodes.each do |start_node|
      log.debug "Trying to assemble from #{start_node.node_id}" if log.debug?

      # If we've already covered this node, don't try it again
      if seen_nodes.include?([start_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]) or
        seen_nodes.include?([start_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST])
        log.debug "Already seen this node, not inspecting further" if log.debug?
        next
      end

      # first attempt to go forward as far as possible, then reverse the path
      # and continue until cannot go farther
      reversed_path_forward = find_beginnning_trail_from_node(start_node, seen_nodes)
      if reversed_path_forward.nil?
        log.debug "Could not find forward path from this node, giving up" if log.debug?
        next
      end
      # Have we already seen this path before?
      #TODO: add in recoherence logic here
      if seen_last_in_path?(reversed_path_forward, seen_nodes)
        log.debug "Already seen the last node of the reversed path forward: #{reversed_path_forward.trail[-1].to_shorthand}, giving up" if log.debug?
        next
      end
      # Assemble ahead again
      log.debug "reversed_path_forward: #{reversed_path_forward.to_shorthand}" if log.debug?
      path, just_visited_onodes = assemble_from(reversed_path_forward)

      # Remove nodes that have already been seen to prevent duplication
      log.debug "Before removing already seen nodes the second time, path was #{path.length} nodes long" if log.debug?
      remove_seen_nodes_from_end_of_path(path, seen_nodes)
      log.debug "After removing already seen nodes the second time, path was #{path.length} nodes long" if log.debug?

      # Add the now seen nodes to the list
      just_visited_onodes.each do |onode_settable|
        seen_nodes << onode_settable
      end

      # Record which nodes have already been visited, so they aren't visited again
      seen_nodes.merge just_visited_onodes
      if @assembly_options[:min_coverage_of_start_nodes]
        # TODO: this could be better by progress += (starting_nodes_just_visited.length)
        progress.increment
      else
        progress.progress += just_visited_onodes.length unless progress.nil?
      end

      if path.length_in_bp < @assembly_options[:min_contig_size]
        log.debug "Path length (#{path.length_in_bp}) less than min_contig_size (#{@assembly_options[:min_contig_size] }), not recording it" if log.debug?
        next
      end
      log.debug "Found a seemingly legitimate path #{path.to_shorthand}" if log.debug?
      if block_given?
        yield path
      else
        paths.push path
      end
    end
    progress.finish unless progress.nil?

    return paths
  end

  def seen_last_in_path?(path, seen_nodes)
    seen_nodes.include?(path[-1].to_settable)
  end

  def gather_starting_nodes
    if @assembly_options[:min_coverage_of_start_nodes] or @assembly_options[:min_length_of_start_nodes]
      starting_nodes = []
      graph.nodes.each do |node|
        if (@assembly_options[:min_coverage_of_start_nodes].nil? or
          node.coverage >= @assembly_options[:min_coverage_of_start_nodes]) and
          (@assembly_options[:min_length_of_start_nodes].nil? or
          node.length_alone >= @assembly_options[:min_length_of_start_nodes])

          starting_nodes.push node
        end
      end
      return starting_nodes
    else
      return graph.nodes
    end
  end

  def setup_progressbar(num_nodes)
    progress = nil
    if @assembly_options[:progressbar_io]
      progress = ProgressBar.create(
        :title => "Assembly",
        :format => '%a %bᗧ%i %p%% %E %t',
        :progress_mark  => ' ',
        :remainder_mark => '･',
        :total => num_nodes,
        :output => @assembly_options[:progressbar_io]
      )
    end
    return progress
  end

  # Given a node, return a path that does not include any short tips, or nil if none is
  # connected to this node.
  # With this path, you can explore forwards. This isn't very clear commenting, but
  # I'm just making this stuff up
  def find_beginnning_trail_from_node(node, previously_seen_nodes)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST #go backwards first, because the path will later be reversed
    dummy_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    dummy_trail.trail = [onode]

    find_node_from_non_short_tip = lambda do |dummy_trail|
      # go all the way forwards
      path, visited_nodes = assemble_from(dummy_trail)

      # Remove already seen nodes from the end of the trail, because
      # they are already included in other paths and this shows
      # up as duplicated contig stretches and this is not correct
      log.debug "Before removing already seen nodes the first time, path was #{path.length} nodes long" if log.debug?
      remove_seen_nodes_from_end_of_path(path, previously_seen_nodes)
      log.debug "After removing already seen nodes the first time, path was #{path.length} nodes long" if log.debug?

      # reverse the path
      path.reverse!
      # peel back up we aren't in a short tip (these lost nodes might be
      # re-added later on)
      cannot_remove_any_more_nodes = false
      log.debug "Before pruning back, trail is #{path.to_shorthand}" if log.debug?
      is_tip, whatever = is_short_tip?(path[-1])
      while is_tip
        if path.length == 1
          cannot_remove_any_more_nodes = true
          break
        end
        path.delete_at(path.length-1)
        log.debug "After pruning back, trail is now #{path.to_shorthand}" if log.debug?
        is_tip, whatever = is_short_tip?(path[-1])
      end

      if cannot_remove_any_more_nodes
        nil
      else
        path
      end
    end

    log.debug "Finding nearest find_connected_node_on_a_path #{node.node_id}" if log.debug?
    if !is_short_tip?(onode)[0]
      log.debug "fwd direction not a short tip, going with that" if log.debug?
      path = find_node_from_non_short_tip.call(dummy_trail)
      if !path.nil?
        return path
      end
    end

    log.debug "rev direction is short tip, now testing reverse" if log.debug?
    onode.reverse!
    if is_short_tip?(onode)[0]
      log.debug "short tip in both directions, there is no good neighbour" if log.debug?
      #short tip in both directions, so not a real contig
      return nil
    else
      log.debug "reverse direction not a short tip, going with that" if log.debug?
      return find_node_from_non_short_tip.call(dummy_trail)
    end
  end

  def remove_seen_nodes_from_end_of_path(path, seen_nodes)
    log.debug "Removing from the end of the path #{path.to_shorthand} any nodes in set of size #{seen_nodes.length}" if log.debug?
    while !path.trail.empty?
      last_node_index = path.length-1
      last_node = path[last_node_index]

      if seen_nodes.include?([last_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]) or
        seen_nodes.include?([last_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST])
        path.trail.delete_at(last_node_index)
      else
        # Last node is not previously seen, chop no further.
        break
      end
    end
    return path
  end

  # Assemble considering reads all reads as single ended. Options:
  # :max_tip_length: if a path is shorter than this in bp, then it will be clipped from the path. Default 100
  # :recoherence_kmer: attempt to separate paths by going back to the reads with this larger kmer
  # :leash_length: don't continue assembly from nodes farther than this distance (in bp) away
  def assemble_from(initial_path, visited_onodes=Set.new)
    options = @assembly_options

    recoherencer = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

    path = initial_path.copy
    visited_onodes = Set.new
    initial_path[0...-1].each do |onode| #Add all except the last node to already seen nodes list
      visited_onodes << onode.to_settable
    end

    dummy_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    oneighbours = nil
    while true
      log.debug "Now assembling from #{path[-1].to_shorthand}" if log.debug?
      if visited_onodes.include?(path[-1].to_settable)
        log.debug "Found circularisation in path, going no further" if log.debug?
        break
      else
        visited_onodes << path[-1].to_settable
      end

      if options[:leash_length] and path.length_in_bp-@graph.hash_length > options[:leash_length]
        log.debug "Beyond leash length, going to further with assembly" if log.debug?
        break
      end

      oneighbours = path.neighbours_of_last_node(@graph)
      if oneighbours.length == 0
        log.debug "Found a dead end, last node is #{path[-1].to_shorthand}" if log.debug?
        break

      elsif oneighbours.length == 1
        to_add = oneighbours[0]
        log.debug "Only one way to go, so going there, to #{to_add.to_shorthand}" if log.debug?
        path.add_oriented_node to_add

      else
        # Reached a fork (or 3 or 4-fork), which way to go?

        # Remove neighbours that are short tips
        oneighbours.reject! do |oneigh|
          is_tip, visiteds = is_short_tip?(oneigh)
          if is_tip
            visiteds.each do |onode_settable|
              visited_onodes << onode_settable
            end
          end
          is_tip
        end

        if oneighbours.length == 0
          log.debug "Found a dead end at a fork, last node is #{path[-1].to_shorthand}" if log.debug?
          break
        elsif oneighbours.length == 1
          log.debug "Clipped short tip(s) off, and then there was only one way to go" if log.debug?
          path.add_oriented_node oneighbours[0]
        elsif options[:recoherence_kmer].nil?
          if log.debug?
            neighbours_string = oneighbours.collect do |oneigh|
              oneigh.to_shorthand
            end.join(' or ')
            log.debug "Came across what appears to be a legitimate fork to nodes #{neighbours_string} and no recoherence kmer given, so giving up" if log.debug?
          end
          break
        else
          unless options[:recoherence_kmer].nil?
            log.debug "Attempting to resolve fork by recoherence" if log.debug?
            oneighbours.select! do |oneigh|
              dummy_trail.trail = path.trail+[oneigh]
              recoherencer.validate_last_node_of_path_by_recoherence(
                dummy_trail,
                options[:recoherence_kmer],
                options[:sequences],
                options[:min_confirming_recoherence_kmer_reads],
                )
            end
          end
          if oneighbours.length == 0
            log.debug "no neighbours passed recoherence, giving up" if log.debug?
            break
          elsif oneighbours.length == 1
            log.debug "After recoherence there's only one way to go, going there"
            path.add_oriented_node oneighbours[0]
          elsif options[:max_coverage_at_fork]
            oneighbours.select! do |oneigh|
              oneigh.node.coverage <= options[:max_coverage_at_fork]
            end
            log.debug "Found #{oneighbours.length} neighbours after removing nodes over max coverage" if log.debug?

            if oneighbours.length == 1
              log.debug "After removing too much coverage neighbours there's only one way to go, going there"
              path.add_oriented_node oneighbours[0]
            else
              log.debug "After removing max coverage nodes, #{oneighbours.length} neighbours found (#{oneighbours.collect{|o| o.to_shorthand}.join(",") }), giving up" if log.debug?
              break
            end


          else
            log.debug "Still forked after recoherence (to #{oneighbours.collect{|on| on.to_shorthand}.join(' & ') }), so seems to be a legitimate fork, giving up" if log.debug?
            break
          end
        end
      end
    end

    visited_onodes << path[-1].to_settable

    return path, visited_onodes
  end

  # Returns false iff there is a path longer than max_tip_length
  # starting at the given oriented_node. Currently works like Dijkstra's
  # shortest path algorithm except that it finds the longest path, not the
  # shortest.
  def is_short_tip?(oriented_node)
    max_tip_length = @assembly_options[:max_tip_length]

    stack = DS::Stack.new
    first = MaxDistancedOrientedNode.new
    first.onode = oriented_node
    first.distance = oriented_node.node.length_alone
    stack.push first

    cache = {}

    while current_max_distanced_onode = stack.pop
      return false, [] if current_max_distanced_onode.distance > max_tip_length

      current_max_distanced_onode.onode.next_neighbours(@graph).each do |oneigh|
        neighbour_distance = current_max_distanced_onode.distance + oneigh.node.length_alone
        next if cache[oneigh.to_settable] and cache[oneigh.to_settable] >= neighbour_distance
        distanced_node = MaxDistancedOrientedNode.new
        distanced_node.onode = oneigh
        distanced_node.distance = neighbour_distance
        cache[oneigh.to_settable] = neighbour_distance
        stack.push distanced_node
      end
    end

    return true, cache.collect{|donode| donode[0]}
  end

  class MaxDistancedOrientedNode
    attr_accessor :onode, :distance
  end
end
