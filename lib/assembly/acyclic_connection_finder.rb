require 'ds'
require 'set'

module Bio
  module AssemblyGraphAlgorithms
    class AcyclicConnectionFinder
      def log
        Bio::Log::LoggerPlus['finishm']
      end

      def find_trails_between_nodes_depth_first_search(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        found_a_path = false
        discovered_list = Set.new

        stack = DS::Stack.new
        stack.push [initial_node, start_looking_off_the_end_of_the_first_node]
        while !stack.empty?
          current_node, current_direction = stack.pop
          if current_node == terminal_node
            log.debug "Found a path between initial and terminal nodes."
            found_a_path = true
          end
        end
        raise
      end

      # Perform a search of the graph starting at the initial node, and try to
      # get all the paths between the initial and terminal nodes. So that things don't get out
      # of control, two bounds on this problem are given:
      #
      # 1. No cyclical trails are followed
      # 2. The total trail length (in base pairs) is less than the leash_length
      #
      # Return an array of Trail objects
      def find_trails_between_nodes(graph, initial_node, terminal_node, leash_length, start_looking_off_the_end_of_the_first_node)
        successful_trails = []
        search_underneath_nodes = lambda do |trail, new_node|
          trail.add_node new_node

          if new_node == terminal_node
            # We have found a successful trail. Done.
            log.info "Successful trail found: #{trail.collect{|node| node.node_id}.join(',')}" if log and log.debug?
            successful_trails.push trail

          else
            # OK so there is more to go (or we have strayed off the golden paths)

            # Need all the neighbours of the current node.
            neighbours = nil
            if trail.last_added_node_dangling_side == :end
              neighbours = graph.neighbours_off_end(new_node)
            elsif trail.last_added_node_dangling_side == :start
              neighbours = graph.neighbours_into_start(new_node)
            else
              raise "Programming error when creating trails between two nodes"
            end

            if neighbours.empty?
              # Dead end path
              log.debug "Dead end path found: #{trail}" if log and log.debug?
            else
              # There is one or more paths from here
              neighbours.each_with_index do |next_node, i|
                if trail.include?(next_node)
                  log.warn "Cyclical trail found in graph, assuming that cyclical path is a dead end." if log and log.debug?
                elsif trail.nucleotide_length + next_node.length_alone > leash_length
                  # We've strayed an unacceptable length away from the thing, forget it
                  log.debug "Straying too far trying to find paths, giving up on this trail. Length with next node is #{trail.nucleotide_length + next_node.length_alone}" if log and log.debug?
                else
                  # Keep on truckin'
                  log.debug "Truckin' onto the next node, from #{new_node.node_id} to #{next_node.node_id}" if log and log.debug?
                  arcs = graph.get_arcs_by_node new_node, next_node
                  arcs.reject!{|arc| arc.begin_node_id == arc.end_node_id}
                  if arcs.length != 1
                    raise "Programming error"
                  end
                  arc = arcs[0]
                  if arc.connects_to_beginning?(next_node.node_id)
                    trail.last_added_node_dangling_side = :end
                  elsif arc.connects_to_end?(next_node.node_id)
                    trail.last_added_node_dangling_side = :start
                  else
                    raise "Programming error"
                  end
                  search_underneath_nodes.call(trail.copy, next_node)
                end
              end
            end
          end
          log.debug "Ending recursion" if log and log.debug?
        end

        # initialise the search
        current_trail = Trail.new
        if start_looking_off_the_end_of_the_first_node
          current_trail.last_added_node_dangling_side = :end
        else
          current_trail.last_added_node_dangling_side = :start
        end
        search_underneath_nodes.call current_trail, initial_node

        return successful_trails
      end

      class Trail < Array
        attr_accessor :last_added_node_dangling_side

        def add_node(node)
          push node
        end

        # Number of nucleotides in this trail
        def nucleotide_length
          return 0 if empty?
          len = self[0].length
          each_with_index do |node, i|
            next if i == 0
            len += node.length_alone
          end
          return len
        end

        # Return a copy of the current trail such that modifying the path
        # of the new node does not affect the path of the original
        def copy
          new_trail = Trail.new
          new_trail.last_added_node_dangling_side = @last_added_node_dangling_side
          each do |node|
            new_trail.push node
          end
          return new_trail
        end

        def to_s
          "Trail #{self.object_id}: #{collect{|n| n.node_id}.join(',')}"
        end
      end
    end
  end
end
