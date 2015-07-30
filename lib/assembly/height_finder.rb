require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::HeightFinder
  include Bio::FinishM::Logging

  # visit nodes in range and determine heights
  def traverse(graph, initial_nodes, options={})
    by_height = []
    traversal_nodes = {}
    cycles = {}
    nodes_in_retrace_phase = Set.new

    # depth-first so stack
    stack = DS::Stack.new
    initial_nodes.each do |onode|
      next if options[:range] and options[:range].none?{|other| other == onode.node }
      traversal_node = CyclicTraversalNode.new
      traversal_node.onode = options[:reverse] ? onode.reverse : onode
      traversal_node.nodes_in = []
      traversal_nodes[traversal_node.onode.to_settable] = traversal_node
      stack.push traversal_node
    end

    while traversal_node = stack.pop
      settable = traversal_node.onode.to_settable
      describe = nil

      if log.debug?
        log.debug "visiting #{traversal_node.describe}."
      end

      # Consider node solved if height is known.
      if not traversal_node.height.nil?
        log.debug "Height of #{traversal_node.describe} is known. Skip." if log.debug?
        next
      end

      # find neighbours
      neighbours = traversal_node.nodes_out
      if neighbours.nil?
        neighbours = traversal_node.onode.next_neighbours(graph)
        if options[:range]
          neighbours.reject!{|onode| options[:range].none?{|other| other == onode.node}} #not in defined range
        end

        # Get or create traversal version of node
        neighbours = neighbours.collect do |onode|
          nbr_settable = onode.to_settable
          traversal_nbr = traversal_nodes[nbr_settable]
          if traversal_nbr.nil?
            traversal_nbr = CyclicTraversalNode.new
            traversal_nbr.onode = onode
            traversal_nbr.nodes_in = []
            traversal_nodes[nbr_settable] = traversal_nbr
          end
          traversal_nbr
        end

        #remember neighbours
        traversal_node.nodes_out = neighbours
      end


      # Can we solve the node?
      if neighbours.empty? #check for a tip
        log.debug "#{traversal_node.describe} is a tip." if log.debug?
        traversal_node.height = 0
        if by_height[0].nil?
          by_height[0] = [traversal_node]
        else
          by_height[0].push(traversal_node)
        end
        log.debug "Found height '0' for #{traversal_node.describe}." if log.debug?
        next
      end

      if nodes_in_retrace_phase.include? settable
        log.debug "Retracing back to #{traversal_node.describe}." if log.debug?

        # Neighbours should have been explored
        # Are neighbours involved in cycles?
        cyclic_neighbours = neighbours.reject{|node| node.cycles.nil?}
        if not cyclic_neighbours.empty?
          # current node is in a cycle if a neighbour is in an unclosed cycle
          log.debug "Found cyclic neighbours #{cyclic_neighbours.collect{|node| node.describe}.join(',')}." if log.debug?
          cyclic_neighbours.each do |node|
            node.cycles.each do |cycle|
              log.debug "Merging cycle #{cycle.onodes.collect{|onode| onode.to_shorthand}.join(',')}." if log.debug?
              new_cycle = traversal_node.merge_unclosed_cycle cycle.copy
              if not new_cycle.nil? and new_cycle.closed?
                log.debug "Cycle completes at #{traversal_node.describe}."
                new_cycle_key = new_cycle.to_settable
                if cycles.has_key? new_cycle_key
                  log.debug "Already seen this cycle." if log.debug?
                else
                  cycles[new_cycle_key] = new_cycle.onodes
                end
              end
            end
          end
        end

        # Unsolved neighbours imply a closed cyclic path.
        # Are neighbours unsolved?
        solved_neighbours = neighbours.reject{|node| node.height.nil?}
        unless solved_neighbours.empty?
          log.debug "We know the heights of neighbours #{solved_neighbours.collect{|node| node.describe}.join(',')}." if log.debug?
          # Compute height from solved neighbours
          height = solved_neighbours.map{|node| node.height}.max + 1
          log.debug "Found height '#{height}' for #{traversal_node.describe}." if log.debug?
          traversal_node.height = height
          if by_height[height].nil?
            by_height[height] = [traversal_node]
          else
            by_height[height].push(traversal_node)
          end
        end
        # If no solved neighbours, leave unsolved

        # Move out of retrace phase
        nodes_in_retrace_phase.delete settable
        log.debug "Finished retracing #{traversal_node.describe}." if log.debug?
        next
      end

      # Move current node to retrace phase, before checking for retracing neighbours in case is own neighbour
      nodes_in_retrace_phase << settable

      # Look for currently retracing neighbours and initiate cycles
      retracing_neighbours = neighbours.select{|node| nodes_in_retrace_phase.include? node.onode.to_settable}
      if not retracing_neighbours.empty?
        log.debug "Initiating cycles for neighbours #{retracing_neighbours.collect{|node| node.describe}.join(',')} currently retracing." if log.debug?
        # initiate cycles for each retracing neighbour
        retracing_neighbours.each{|node| traversal_node.initiate_cycle(node.onode)}
      end

      # Return node stack and push neighbours
      stack.push traversal_node
      log.debug "Pushing #{traversal_node.describe} in retrace mode." if log.debug?
      neighbours.each do |node|
        node_settable = node.onode.to_settable

        # Note the parent of neighbour unless already known
        nodes_in = node.nodes_in
        if not nodes_in.any?{|nbr| nbr.onode == node.onode}
          nodes_in.push traversal_node
        end

        if nodes_in_retrace_phase.include? node_settable
          # A currently retracing neighbour implies a cycle, cut it off here
          log.debug "Neighbour #{node.describe} is retracing. Not revisiting." if log.debug?
        else
          log.debug "Pushing neighbour #{node.describe}." if log.debug?
          stack.push node
        end
      end
    end
    return by_height, cycles.values
  end

  class TraversalNode
    attr_accessor :onode, :height, :nodes_in, :nodes_out

    def describe
      @onode.to_shorthand
    end

    def node_id
      @onode.node_id
    end
  end

  class CyclicTraversalNode < TraversalNode
    attr_accessor :cycles

    def initiate_cycle(onode)
      cycle = CyclePath.new
      cycle.onodes = [onode]
      merge_unclosed_cycle cycle
    end

    def merge_unclosed_cycle(cycle)
      return if cycle.closed?
      if cycle.onodes.last == @onode
        cycle.closed = true
      else
        cycle.onodes.unshift @onode
      end
      if @cycles.nil?
        @cycles = [cycle]
      else
        @cycles.push(cycle)
      end
      return cycle
    end

    class CyclePath
      attr_accessor :onodes, :closed

      def closed?
        return @closed == true
      end

      def copy
        cycle = CyclePath.new
        cycle.onodes = @onodes[0..-1]
        cycle.closed = @closed
        cycle
      end

      def to_settable
        # return sorted list of onode settables
        @onodes.collect{|onode| onode.to_settable}.sort do |a, b|
          result = a[0] <=> b[0]
          if result == 0
            result = a[1] == Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode::START_IS_FIRST ? -1 : 1
          end
          result
        end
      end
    end
  end


  # maximum paths
  def max_paths_through(by_height)
    max_paths_from = {}
    by_height.each_with_index do |nodes, level|
      log.debug "At height #{level}." if log.debug?
      if level == 0  # tips
        nodes.each do |node|
          log.debug "Counted maximum of 1 path to #{node.describe}." if log.debug?
          max_paths_from[node.onode.to_settable] = 1
        end
        next
      end

      nodes.each do |node|
        settable = node.onode.to_settable
        max_paths_from_neighbours = node.nodes_out.collect{|nbr| max_paths_from[nbr.onode.to_settable]}.reject{|n| n.nil?}
        log.debug "Found neighbours of #{node.describe} with maximum paths #{max_paths_from_neighbours.join(',')}." if log.debug?
        max_paths_from[settable] = max_paths_from_neighbours.reduce{|memo, num| memo+num}
        log.debug "Counted maximum of #{max_paths_from[settable]} paths to #{node.describe}." if log.debug?
      end
    end

    # Get the graph roots (which are nodes with no parents) and add max_paths_from for each to get graph total
    root_keys = by_height.flatten.select{|node| node.nodes_in.empty? }.collect{|node| node.onode.to_settable}
    log.debug "Found graph roots #{root_keys.collect{|settable| settable[0]}.join(',')} with maximum paths #{root_keys.collect{|key| max_paths_from[key]}.join(',')}." if log.debug?
    max_paths = root_keys.map{|settable| max_paths_from[settable]}.reduce{|memo, num| memo+num}
    log.debug "Counted maximum of #{max_paths} through graph." if log.debug?
    return max_paths
  end

  # minimum paths
  def min_paths_through(by_height)
    live_nodes = Set.new
    max_alive_counter = 0
    by_height.each_with_index do |nodes, level|
      log.debug "At height #{level}." if log.debug?
      # nodes at current level become live
      nodes.each do |node|
        settable = node.onode.to_settable
        log.debug "Setting #{node.describe} as live." if log.debug?
        live_nodes << settable
      end
      if level > 0
        #children of nodes at current level are no longer live
        nodes.each do |node|
          children = node.nodes_out
          children.each do |nbr|
            log.debug "Setting child #{nbr.describe} of live node #{node.describe} as inactive." if log.debug?
            live_nodes.delete(nbr.onode.to_settable)
          end
        end
      end

      log.debug "There are currently #{live_nodes.length} nodes alive. Max is #{max_alive_counter}." if log.debug?
      if live_nodes.length > max_alive_counter
        #track the maximum live nodes at any level
        log.debug "Updating max to #{live_nodes.length}." if log.debug?
        max_alive_counter = live_nodes.length
      end
    end
    return max_alive_counter
  end

  def find_oriented_edge_of_range(graph, nodes=nil)
    nodes ||= graph.nodes
    log.debug "Looking for oriented start and end points from #{nodes.collect{|n| n.node_id}.join(',')}" if log.debug?
    nodes_all_directions = nodes.collect{|node| [[node, true], [node, false]]}.flatten(1)


    # Find nodes and directions which are not reachable from other nodes within range
    unreached_nodes = {}
    nodes_all_directions.each do |node_and_direction|
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new node_and_direction[0], node_and_direction[1]
      unless unreached_nodes.has_key? onode.to_settable
        unreached_nodes[onode.to_settable] = onode
      end
      onode.next_neighbours(graph).each do |oneigh|
        unreached_nodes[oneigh.to_settable] = nil
      end
    end

    entry_points = unreached_nodes.values.reject{|n| n.nil?}
    log.debug "Found the following nodes for a particular orientation have no paths connecting to other nodes in range: #{entry_points.collect{|n| n.to_shorthand}.join(',')}" if log.debug?

    # Start from an unreachable node, and trace all paths until the reverse end of other unreachable nodes
    # are reached, which are then defined as 'end' nodes. When finished, choose a remaining non-end unreachable
    # node and repeat, stopping paths if an already seen node is encountered.
    seen_nodes = Set.new
    start_onodes = []
    end_onodes = []
    stack = DS::Stack.new
    entry_points.reverse.each do |onode|
      stack.push onode
    end

    while current_node = stack.pop
      log.debug "At node #{current_node.to_shorthand}" if log.debug?

      node_id = current_node.node_id
      if seen_nodes.include? node_id or not nodes.include? current_node.node
        log.debug "Node has been seen or is out of range. Skipping..." if log.debug?
        next
      end
      seen_nodes << node_id

      current_unreached = unreached_nodes[current_node.to_settable]
      log.debug "Is current unreached? #{current_unreached}" if log.debug?
      if current_unreached
        log.debug "Defining starting node #{current_unreached.to_shorthand}" if log.debug?
        # Found start node
        start_onodes.push current_unreached
      else
        reverse_unreached = unreached_nodes[current_node.reverse.to_settable]
        log.debug "Is reverse unreached? #{reverse_unreached}" if log.debug?
        if reverse_unreached
          log.debug "Found ending node #{reverse_unreached.to_shorthand}" if log.debug?
          # Found end node
          end_onodes.push reverse_unreached
        end
      end

      current_node.next_neighbours(graph).each do |onode|
        log.debug "Adding neighbour #{onode.to_shorthand} to stack" if log.debug?
        stack.push onode
      end
    end

    return start_onodes, end_onodes
  end
end
