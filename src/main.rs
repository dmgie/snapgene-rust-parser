use std::fs::File;

use nom::{
    bytes::complete::{tag, take},
    combinator::{map_res, opt},
    sequence::{pair, preceded,separated_pair},
    // bits::{bits, bytes, streaming::take},
    error::Error,
    IResult,
};


// Snapgene has "packet" types i.e sections
// The structure is as follows:
// 1 byte packet type (determines the rest of the packet structure)
// 4 byte (big-endian long int) packet length (in bytes)
// \inf number of data bytes
// #[derive(Debug, PartialEq)]
// struct SnapgenePacket {
//     packet_tag: PacketType,
//     packet_length: u32,
//     data: Vec<u8>,
// }
type Long = u32;
type Short = u16;

// A packet header is the first 5 bytes of a packet
#[derive(Debug, PartialEq)]
struct PacketHeader {
    packet_tag: PacketType,
    packet_length: u32,
}

impl PacketHeader {
    fn parse(input: &[u8]) -> PacketHeader {
        let parsed = PacketHeader::parse_header(input);
        match parsed {
            Ok((_, snap_parser)) => snap_parser,
            Err(e) => panic!("Error parsing snapgene packet header : {:?}", e),
        }
    }

    fn parse_header(input: &[u8]) -> IResult<&[u8], PacketHeader> {
        // Take first byte
        let (input, packet_tag) = take(1usize)(input)?;

        // Take next 4 bytes, long int
        let (input, packet_length) = take(4usize)(input)?;

        // Set the values
        let packet_header = PacketHeader {
            packet_tag: match packet_tag[0] {
                0x09 => PacketType::Cookie,
                0x00 => PacketType::Sequence,
                0x0A => PacketType::Feature,
                0x06 => PacketType::Notes,
                _ => panic!("Unknown packet type"),
            },
            packet_length: u32::from_be_bytes(packet_length.try_into().unwrap()),
        };

        Ok((input, packet_header))
    }
}

// Cookie is the initial packet
#[derive(Debug, PartialEq)]
struct Cookie {
    packet_tag: PacketType,
    packet_length: u32,
    magic_cookie: [u8; 8],
    seq_type: SeqType,
    export_version: u16,
    import_version: u16,
}

// Only appears once
// Sequence length = packet length - 1 since the 1st byte for the sequence_flag
#[derive(Debug, PartialEq)]
struct SequencePacket {
    packet_tag: PacketType,
    packet_length: u32,
    sequence_flag: Topology,
    seq: Vec<u8>,
}

#[derive(Debug, PartialEq)]
enum Topology {
    Linear = 0x00, // Absent flag = linear
    Circular, // Otherwise circular
}

#[derive(Debug, PartialEq)]
enum PacketType {
    Cookie = 0x09,
    Sequence = 0x00,
    Notes = 0x06,
    Feature = 0x0A,
}

#[derive(Debug, PartialEq)]
enum SeqType {
    DNA = 1,
    RNA = 2,
    Protein = 3,
}

impl SequencePacket {
    fn parse(input: &[u8]) -> SequencePacket {
        let parsed = SequencePacket::parse_dnapacket(input);
        match parsed {
            Ok((_, snap_parser)) => snap_parser,
            Err(e) => panic!("Error parsing snapgene file: {:?}", e),
        }
    }

    fn parse_dnapacket(input: &[u8]) -> IResult<&[u8], SequencePacket> {
        // Take first byte
        let (input, packet_tag) = take(1usize)(input)?;

        // Take next 4 bytes i.e a long
        let (input, packet_length) = take(4usize)(input)?;

        // Take next byte
        let (input, sequence_flag) = take(1usize)(input)?;

        // Get actual sequence: packet_length - 1 (sequence_flag)
        let length = u32::from_be_bytes(packet_length.try_into().unwrap()) - 1;
        let (input, seq) = take(length as usize)(input)?;



        // Set the values
        let snap_parser = SequencePacket {
            packet_tag: match packet_tag[0] {
                0x09 => PacketType::Cookie,
                0x00 => PacketType::Sequence,
                0x0A => PacketType::Feature,
                0x06 => PacketType::Notes,
                _ => panic!("Unknown packet type"),
            },
            packet_length: u32::from_be_bytes(packet_length.try_into().unwrap()),
            sequence_flag: match sequence_flag[0] {
                0x00 => Topology::Linear,
                _ => Topology::Circular,
            },
            seq: seq.to_vec(),
        };

        Ok((input, snap_parser))
    }
}

// This is a binary parser
impl Cookie {
    fn parse(input: &[u8]) -> Cookie {
        let parsed = Cookie::parse_cookie(input);
        match parsed {
            Ok((_, snap_parser)) => snap_parser,
            Err(e) => panic!("Error parsing snapgene file: {:?}", e),
        }
    }
    // Cookie always starts with 0x09
    fn parse_cookie(input: &[u8]) -> IResult<&[u8], Cookie> {
        // Take first byte
        let (input, packet_tag) = take(1usize)(input)?;

        // Take next 4 bytes, long int
        let (input, packet_length) = take(4usize)(input)?;

        // Take next 8 bytes, string
        let (input, magic_cookie) = take(8usize)(input)?;

        // Take next 2 bytes, a short integer
        let (input, seq_type) = take(2usize)(input)?;

        let (input, export_version) = take(2usize)(input)?;

        let (input, import_version) = take(2usize)(input)?;

        // Set the values
        let cookie_packet = Cookie {
            packet_tag: match packet_tag[0] {
                0x09 => PacketType::Cookie,
                0x00 => PacketType::Sequence,
                0x0A => PacketType::Feature,
                0x06 => PacketType::Notes,
                _ => panic!("Unknown packet type"),
            },
            packet_length: u32::from_be_bytes(packet_length.try_into().unwrap()),
            magic_cookie: magic_cookie.try_into().unwrap(),
            seq_type: match u16::from_be_bytes(seq_type.try_into().unwrap()) {
                1 => SeqType::DNA,
                2 => SeqType::RNA,
                3 => SeqType::Protein,
                _ => panic!("Unknown sequence type"),
            },
            export_version: u16::from_be_bytes(export_version.try_into().unwrap()),
            import_version: u16::from_be_bytes(import_version.try_into().unwrap()),
        };

        Ok((input, cookie_packet))
    }
}

fn main() {
    let f = std::fs::read("./pcDNA4 TO.dna").unwrap();
    let (input, cookie_packet) = Cookie::parse_cookie(&f).unwrap();
    println!("{cookie_packet:?}");
    println!("{:?}", &input[0..15]);
    let (input, dna) = SequencePacket::parse_dnapacket(input).unwrap();
    println!("{dna:?}");
    let dna = String::from_utf8(dna.seq).unwrap();

    // Take the next 4 bytes, long int
    let (input, next_packet) = take_long(input).unwrap();
    println!("{next_packet:?}");

}
